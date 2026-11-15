#!/usr/bin/env python3
"""Independent diagnosis for the pointwise trend batch.

All expectations in this file are rebuilt from the raw observations in
data/*.json with exact rational arithmetic.  The original processing
chain (resample.materialize, summary.summarize/merge, trend.fit) is
NOT used to derive any expected value; it is only executed as the
"original entry" whose output is compared against the independent model
(and whose two incident failures are recorded).

Only the neutral exact-arithmetic helpers in matrix.py are reused
(rank / inverse / multiply / transpose), as CONTRACT.md permits.

Modes
  build   regenerate evidence.json (not part of the run manifest)
  verify  recompute all evidence, check REPORT.md anchors, run the
          original normal/local regressions, run the original incident
          entry twice (both must fail with TREND_PARTITION_DRIFT), and
          cross-check every plan against the original pipeline entry
  demo    show the canonical batch fit, the existing fit and a genuine
          single-factor contrast

Neither verify nor demo writes to data/, logs/ or any original file.
"""

import argparse
import copy
import hashlib
import json
import math
import subprocess
import sys
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from matrix import rank, inverse, multiply, transpose  # neutral arithmetic only

MARKER = "TREND_PARTITION_DRIFT"
PLANS = ["all", "cut", "pairs", "striped", "single"]
TOL = 1e-10


# --------------------------------------------------------------------------
# exact parsing / encoding
# --------------------------------------------------------------------------
def exact(q):
    """Parse the contract's integer or [numerator, positive denominator]."""
    if isinstance(q, list):
        if len(q) != 2 or not isinstance(q[1], int) or isinstance(q[1], bool) or q[1] <= 0:
            raise ValueError("bad rational: %r" % (q,))
        return F(q[0], q[1])
    if isinstance(q, bool) or not isinstance(q, int):
        raise ValueError("numbers must be integers or [numerator, positive denominator]")
    return F(q)


def frac(x):
    return [x.numerator, x.denominator]


def encode(obj):
    if isinstance(obj, F):
        return frac(obj)
    if isinstance(obj, dict):
        return {k: encode(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [encode(v) for v in obj]
    return obj


def read_json(path):
    return json.loads(Path(path).read_text())


def parse_package(pkg):
    """Contract-semantics materialisation, written independently."""
    source = {}
    for row in pkg["observations"]:
        ident = row["id"]
        if ident in source:
            raise ValueError("duplicate source id: %s" % ident)
        item = {k: exact(row[k]) for k in ("time", "value", "variance")}
        if item["variance"] <= 0:
            raise ValueError("nonpositive variance: %s" % ident)
        source[ident] = item
    windows = []
    seen = set()
    for spec in pkg["windows"]:
        if spec["id"] in seen:
            raise ValueError("duplicate window id: %s" % spec["id"])
        seen.add(spec["id"])
        weights = {k: exact(v) for k, v in spec["weights"].items()}
        if not weights or any(k not in source or v < 0 for k, v in weights.items()):
            raise ValueError("invalid weight in %s" % spec["id"])
        if sum(weights.values()) != 1:
            raise ValueError("weights of %s must sum to one" % spec["id"])
        weights = {k: v for k, v in weights.items() if v}
        windows.append({
            "id": spec["id"],
            "weights": weights,
            "time": sum(v * source[k]["time"] for k, v in weights.items()),
            "value": sum(v * source[k]["value"] for k, v in weights.items()),
        })
    plans = {name: [[wid for wid in group] for group in blocks]
             for name, blocks in pkg["plans"].items()}
    return source, windows, plans


# --------------------------------------------------------------------------
# exact batch model
# --------------------------------------------------------------------------
class Batch:
    """Window weight matrix and known source covariance."""

    def __init__(self, pkg, label):
        self.label = label
        self.source, self.windows, self.plans = parse_package(pkg)
        self.ids = sorted(self.source)
        self.n = len(self.ids)
        self.v = [self.source[k]["variance"] for k in self.ids]
        self.by_id = {w["id"]: w for w in self.windows}
        self.order = [w["id"] for w in self.windows]
        self.rows = {w["id"]: [w["weights"].get(k, F(0)) for k in self.ids]
                     for w in self.windows}

    # window covariance = W diag(v) W^T
    def window_cov(self, a, b):
        wa, wb = self.rows[a], self.rows[b]
        return sum((x * y * v for x, y, v in zip(wa, wb, self.v)), F(0))

    def full_cov(self, wids=None):
        wids = self.order if wids is None else wids
        return [[self.window_cov(a, b) for b in wids] for a in wids]

    def window_variance(self, wid):
        return sum((x * x * v for x, v in zip(self.rows[wid], self.v)), F(0))

    def greedy_basis(self, wids):
        """Block-local selection, the rule used by summary.summarize."""
        rows, selected = [], []
        for wid in wids:
            row = self.rows[wid]
            if rank(rows + [row]) > len(rows):
                rows.append(row)
                selected.append(wid)
        return selected

    def maximal_bases(self):
        """All maximal linearly independent subsets of all window rows."""
        m = len(self.order)
        for size in range(m, 0, -1):
            found = []
            for comb in combinations(self.order, size):
                if rank([self.rows[w] for w in comb]) == size:
                    found.append(list(comb))
            if found:
                return found, size
        return [], 0

    def gls(self, basis):
        """Exact GLS on the given independent window basis."""
        cov = self.full_cov(basis)
        precision = inverse(cov)
        design = [[F(1), self.by_id[w]["time"]] for w in basis]
        observed = [[self.by_id[w]["value"]] for w in basis]
        score = multiply(transpose(design), precision)
        information = multiply(score, design)
        rhs = multiply(score, observed)
        information_inv = inverse(information)
        beta = multiply(information_inv, rhs)
        # source-loading matrix: beta_hat = B * y_source
        basis_rows = [self.rows[w] for w in basis]
        loadings = multiply(information_inv, multiply(score, basis_rows))
        return {
            "basis": basis,
            "covariance_of_basis": cov,
            "precision": precision,
            "design": design,
            "information": information,
            "rhs": rhs,
            "covariance": information_inv,
            "beta": beta,
            "loadings": loadings,
        }

    def source_loadings(self, loadings):
        return {
            "intercept": {k: loadings[0][i] for i, k in enumerate(self.ids)},
            "slope": {k: loadings[1][i] for i, k in enumerate(self.ids)},
        }

    def actual_variance(self, loadings, which):
        row = loadings[0] if which == "intercept" else loadings[1]
        return sum((v * x * x for v, x in zip(self.v, row)), F(0))

    def nullspace(self, wids):
        """A basis of the cross-window linear dependencies (exact)."""
        rows = [self.rows[w] for w in wids]
        m, n = len(rows), len(rows[0])
        a = [[rows[i][j] for i in range(m)] for j in range(n)]  # n x m
        pivots, r = [], 0
        for c in range(m):
            p = next((i for i in range(r, n) if a[i][c]), None)
            if p is None:
                continue
            a[r], a[p] = a[p], a[r]
            d = a[r][c]
            a[r] = [x / d for x in a[r]]
            for i in range(n):
                if i != r:
                    d = a[i][c]
                    a[i] = [x - d * y for x, y in zip(a[i], a[r])]
            pivots.append(c)
            r += 1
            if r == n:
                break
        vectors = []
        for free in (c for c in range(m) if c not in pivots):
            x = [F(0)] * m
            x[free] = F(1)
            for i, c in enumerate(pivots):
                x[c] = -a[i][free]
            vectors.append({wids[i]: frac(x[i]) for i in range(m) if x[i]})
        return vectors

    def canonical(self):
        bases, size = self.maximal_bases()
        chosen = self.gls(bases[0])
        # every maximal independent basis must give the same BLUE
        equivalent = True
        for basis in bases[1:]:
            other = self.gls(basis)
            same = (other["beta"] == chosen["beta"]
                    and other["covariance"] == chosen["covariance"]
                    and other["loadings"] == chosen["loadings"])
            equivalent = equivalent and same
        cov = chosen["covariance"]
        var_i = self.actual_variance(chosen["loadings"], "intercept")
        var_s = self.actual_variance(chosen["loadings"], "slope")
        checks = {
            "intercept_variance_by_loadings": frac(var_i),
            "slope_variance_by_loadings": frac(var_s),
            "covariance_matches_loadings": (var_i == cov[0][0] and var_s == cov[1][1]),
            "sum_slope_loadings_zero": sum(chosen["loadings"][1]) == 0,
            "sum_slope_loadings_time_one": sum(
                (x * self.source[k]["time"]
                 for x, k in zip(chosen["loadings"][1], self.ids)), F(0)) == 1,
        }
        result = self.report(
            chosen["beta"], cov[1][1], size, len(self.windows))
        return {
            "basis": chosen["basis"],
            "basis_count": len(bases),
            "maximal_independent_size": size,
            "all_bases_equivalent": equivalent,
            "full_covariance_rank": rank(self.full_cov()),
            "full_covariance": encode(self.full_cov()),
            "nullspace_of_window_rows": self.nullspace(self.order),
            "window_variance": {w: frac(self.window_variance(w)) for w in self.order},
            "window_time": {w: frac(self.by_id[w]["time"]) for w in self.order},
            "window_value": {w: frac(self.by_id[w]["value"]) for w in self.order},
            "information": chosen["information"],
            "rhs": chosen["rhs"],
            "covariance_of_basis": chosen["covariance_of_basis"],
            "covariance": chosen["covariance"],
            "beta": chosen["beta"],
            "loadings": self.source_loadings(chosen["loadings"]),
            "loadings_matrix": chosen["loadings"],
            "variance_checks": checks,
            "result": result,
        }

    def report(self, beta, slope_var, information_rank, window_count):
        se = math.sqrt(float(slope_var))
        slope = beta[1][0]
        return {
            "intercept": beta[0][0],
            "slope": slope,
            "slope_variance": slope_var,
            "standard_error": se,
            "ci95": [float(slope) - 1.96 * se, float(slope) + 1.96 * se],
            "information_rank": information_rank,
            "window_count": window_count,
        }
    # ----------------------------------------------------------------------
    # model of the estimator actually implemented by pipeline.py
    # ----------------------------------------------------------------------
    def existing_blocks(self, plan):
        blocks = []
        for group in plan:
            basis = self.greedy_basis(group)
            cov = self.full_cov(basis)
            precision = inverse(cov)
            design = [[F(1), self.by_id[w]["time"]] for w in basis]
            observed = [[self.by_id[w]["value"]] for w in basis]
            score = multiply(transpose(design), precision)
            blocks.append({
                "group": list(group),
                "basis": basis,
                "block_covariance": cov,
                "precision": precision,
                "design": design,
                "information": multiply(score, design),
                "rhs": multiply(score, observed),
                "loadings": multiply(score, [self.rows[w] for w in basis]),
            })
        return blocks

    def cross_block(self, blocks):
        """Nonzero covariances between basis windows placed in different blocks."""
        owner = {}
        for index, block in enumerate(blocks):
            for w in block["basis"]:
                owner[w] = index
        pairs = []
        for a in self.order:
            for b in self.order:
                if a in owner and b in owner and owner[a] < owner[b]:
                    cov = self.window_cov(a, b)
                    if cov:
                        pairs.append({"windows": [a, b], "covariance": frac(cov)})
        return pairs

    def existing(self, plan):
        blocks = self.existing_blocks(plan)
        info = [[sum((b["information"][i][j] for b in blocks), F(0))
                 for j in range(2)] for i in range(2)]
        rhs = [[sum((b["rhs"][i][0] for b in blocks), F(0))] for i in range(2)]
        info_inv = inverse(info)
        beta = multiply(info_inv, rhs)
        # beta_hat = inv(info) * sum_b S_b W_b * y_source
        total_load = [[sum((b["loadings"][i][k] for b in blocks), F(0))
                       for k in range(self.n)] for i in range(2)]
        coeff = multiply(info_inv, total_load)
        actual_i = self.actual_variance(coeff, "intercept")
        actual_s = self.actual_variance(coeff, "slope")
        reported_var = info_inv[1][1]
        reported = self.report(beta, reported_var,
                               sum(len(b["basis"]) for b in blocks),
                               len(self.windows))
        actual_se = math.sqrt(float(actual_s))
        slope = beta[1][0]

        # Independent cross-check of the actual variance: express the merged
        # estimator as a linear form in the *window* values, then use the full
        # window covariance (which knows about cross-block dependence too).
        basis_windows = [w for b in blocks for w in b["basis"]]
        window_coeff = [[], []]
        for b in blocks:
            local = multiply(info_inv, multiply(
                transpose(b["design"]), b["precision"]))
            window_coeff[0].extend(local[0])
            window_coeff[1].extend(local[1])
        full = self.full_cov(basis_windows)
        cov_sum = sum((window_coeff[0][i] * full[i][j] * window_coeff[0][j]
                       for i in range(len(basis_windows))
                       for j in range(len(basis_windows))), F(0))
        slope_sum = sum((window_coeff[1][i] * full[i][j] * window_coeff[1][j]
                         for i in range(len(basis_windows))
                         for j in range(len(basis_windows))), F(0))

        return {
            "blocks": blocks,
            "information": info,
            "rhs": rhs,
            "information_inverse": info_inv,
            "beta": beta,
            "slope_loadings": {k: coeff[1][i] for i, k in enumerate(self.ids)},
            "intercept_loadings": {k: coeff[0][i] for i, k in enumerate(self.ids)},
            "loadings_matrix": coeff,
            "actual_intercept_variance": actual_i,
            "actual_slope_variance": actual_s,
            "variance_ratio": actual_s / reported_var,
            "actual_standard_error": actual_se,
            "actual_ci95": [float(slope) - 1.96 * actual_se,
                            float(slope) + 1.96 * actual_se],
            "window_covariance_check": {
                "intercept_variance": frac(cov_sum),
                "slope_variance": frac(slope_sum),
                "slope_matches_source_loadings": slope_sum == actual_s,
                "intercept_matches_source_loadings": cov_sum == actual_i,
            },
            "union_window_rank": rank([self.rows[w] for w in basis_windows]),
            "union_nullspace": self.nullspace(basis_windows),
            "cross_block_covariance": self.cross_block(blocks),
            "reported_information_equals_marginal_sum": (
                info == self.marginal_information_sum()),
            "marginal_information_sum": encode(self.marginal_information_sum()),
            "reported": reported,
        }

    def marginal_information_sum(self):
        """Sum over all windows of x x^T / var(w): the 'all independent' model."""
        total = [[F(0), F(0)], [F(0), F(0)]]
        for w in self.order:
            x = [F(1), self.by_id[w]["time"]]
            var = self.window_variance(w)
            for i in range(2):
                for j in range(2):
                    total[i][j] += x[i] * x[j] / var
        return total


# --------------------------------------------------------------------------
# evidence assembly
# --------------------------------------------------------------------------
def decomposition_evidence(batch, canon):
    """Exact split of the merged information/rhs into canonical + phantom gain."""
    jc, rc, bc = canon["information"], canon["rhs"], canon["beta"]
    out = {}
    for name in PLANS:
        ex = batch.existing(batch.plans[name])
        dj = [[ex["information"][i][j] - jc[i][j] for j in range(2)]
              for i in range(2)]
        dr = [[ex["rhs"][i][0] - rc[i][0]] for i in range(2)]
        # along-line score and the orthogonal residual that actually moves beta
        along = [[sum(dj[i][j] * bc[j][0] for j in range(2))] for i in range(2)]
        residual = [[dr[i][0] - along[i][0]] for i in range(2)]
        shift = multiply(inverse(ex["information"]), residual)
        # exact PSD test for a symmetric 2x2: a>=0, c>=0, det>=0
        psd = dj[0][0] >= 0 and dj[1][1] >= 0 and (
            dj[0][0] * dj[1][1] - dj[0][1] * dj[1][0]) >= 0
        out[name] = {
            "information_difference": encode(dj),
            "rhs_difference": encode(dr),
            "information_difference_rank": rank(dj),
            "information_difference_trace": frac(dj[0][0] + dj[1][1]),
            "information_difference_determinant": frac(
                dj[0][0] * dj[1][1] - dj[0][1] * dj[1][0]),
            "information_difference_psd": psd,
            "rank_deficit": ex["reported"]["information_rank"]
                            - canon["result"]["information_rank"],
            "score_along_canonical_line": encode(along),
            "score_orthogonal_residual": encode(residual),
            "implied_estimate_shift": encode(shift),
            "actual_estimate_shift": encode(
                [[ex["beta"][i][0] - bc[i][0]] for i in range(2)]),
            "shift_matches": shift == [[ex["beta"][i][0] - bc[i][0]] for i in range(2)],
        }
    # The cut plan's gain is exactly one extra copy of source c's own
    # observation: dj = x x^T, dr = 5 x with x = (1, 2), variance(c) = 1.
    cut = out["cut"]
    dj = [[exact(q) for q in row] for row in cut["information_difference"]]
    dr = [[exact(q) for q in row] for row in cut["rhs_difference"]]
    out["cut"]["phantom_is_source_c_squared"] = (
        dj == [[F(1), F(2)], [F(2), F(4)]]
        and dr == [[F(5)], [F(10)]]
        and batch.source["c"]["variance"] == 1
        and batch.source["c"]["time"] == 2
        and batch.source["c"]["value"] == 5)
    return out


def duplicate_package(pkg):
    dup = copy.deepcopy(pkg)
    dup["windows"].append({"id": "w8", "weights": {"b": [1, 2], "c": [1, 2]}})
    return dup


def independent_package(pkg):
    """Same numbers, but every window references fresh, independent ids."""
    observations, windows = [], []
    for spec in pkg["windows"]:
        weights = {}
        for k, weight in spec["weights"].items():
            origin = next(r for r in pkg["observations"] if r["id"] == k)
            fresh = "%s::%s" % (spec["id"], k)
            observations.append({
                "id": fresh, "time": origin["time"],
                "value": origin["value"], "variance": origin["variance"],
            })
            weights[fresh] = weight
        windows.append({"id": spec["id"], "weights": weights})
    return {"observations": observations, "windows": windows,
            "plans": copy.deepcopy(pkg["plans"])}


def plan_evidence(batch, canonical):
    out = {}
    for name in PLANS:
        ex = batch.existing(batch.plans[name])
        need(ex["window_covariance_check"]["slope_matches_source_loadings"]
             and ex["window_covariance_check"]["intercept_matches_source_loadings"],
             "window-covariance cross-check failed for plan %s" % name)
        out[name] = {
            "result": encode(ex["reported"]),
            "slope": frac(ex["beta"][1][0]),
            "intercept": frac(ex["beta"][0][0]),
            "block_bases": [b["basis"] for b in ex["blocks"]],
            "block_information": [encode(b["information"]) for b in ex["blocks"]],
            "block_rhs": [encode(b["rhs"]) for b in ex["blocks"]],
            "merged_information": encode(ex["information"]),
            "merged_rhs": encode(ex["rhs"]),
            "merged_information_inverse": encode(ex["information_inverse"]),
            "slope_loadings": {k: frac(v) for k, v in ex["slope_loadings"].items()},
            "intercept_loadings": {k: frac(v) for k, v in ex["intercept_loadings"].items()},
            "reported_slope_variance": frac(ex["reported"]["slope_variance"]),
            "actual_slope_variance": frac(ex["actual_slope_variance"]),
            "actual_intercept_variance": frac(ex["actual_intercept_variance"]),
            "variance_ratio": frac(ex["variance_ratio"]),
            "variance_ratio_float": float(ex["variance_ratio"]),
            "actual_standard_error": ex["actual_standard_error"],
            "actual_ci95": ex["actual_ci95"],
            "window_covariance_check": ex["window_covariance_check"],
            "union_window_rank": ex["union_window_rank"],
            "union_nullspace": ex["union_nullspace"],
            "cross_block_covariance": ex["cross_block_covariance"],
            "reported_information_equals_marginal_sum":
                ex["reported_information_equals_marginal_sum"],
            "information_rank": ex["reported"]["information_rank"],
            "window_count": ex["reported"]["window_count"],
            "delta_vs_canonical": {
                "slope_difference": frac(ex["beta"][1][0] - canonical["beta"][1][0]),
                "intercept_difference": frac(ex["beta"][0][0] - canonical["beta"][0][0]),
                "reported_variance_over_canonical": frac(
                    ex["reported"]["slope_variance"] / canonical["result"]["slope_variance"]),
                "actual_variance_over_canonical": frac(
                    ex["actual_slope_variance"] / canonical["result"]["slope_variance"]),
            },
        }
    return out


def canonical_evidence(canon):
    return {
        "basis": canon["basis"],
        "basis_count": canon["basis_count"],
        "maximal_independent_size": canon["maximal_independent_size"],
        "all_bases_equivalent": canon["all_bases_equivalent"],
        "full_covariance_rank": canon["full_covariance_rank"],
        "full_covariance": encode(canon["full_covariance"]),
        "nullspace_of_window_rows": encode(canon["nullspace_of_window_rows"]),
        "window_variance": encode(canon["window_variance"]),
        "window_time": encode(canon["window_time"]),
        "window_value": encode(canon["window_value"]),
        "information": encode(canon["information"]),
        "rhs": encode(canon["rhs"]),
        "covariance_of_basis": encode(canon["covariance_of_basis"]),
        "covariance": encode(canon["covariance"]),
        "intercept": frac(canon["beta"][0][0]),
        "slope": frac(canon["beta"][1][0]),
        "loadings": encode(canon["loadings"]),
        "variance_checks": encode(canon["variance_checks"]),
        "result": encode(canon["result"]),
    }


def run_command(argv, cwd, stdin=None):
    proc = subprocess.run(argv, cwd=str(cwd), capture_output=True, text=True,
                          timeout=120)
    return {
        "argv": list(argv),
        "exit_code": proc.returncode,
        "stdout": proc.stdout,
        "stderr": proc.stderr,
        "stdout_sha256": hashlib.sha256(proc.stdout.encode()).hexdigest(),
        "stderr_sha256": hashlib.sha256(proc.stderr.encode()).hexdigest(),
    }


def regression_evidence():
    normal = run_command([sys.executable, "checks.py", "normal"], ROOT)
    local = run_command([sys.executable, "checks.py", "local"], ROOT)
    return {
        "normal": {
            "exit_code": normal["exit_code"],
            "result": json.loads(normal["stdout"])["normal"],
            "stdout_sha256": normal["stdout_sha256"],
        },
        "local": {
            "exit_code": local["exit_code"],
            "result": json.loads(local["stdout"]),
            "stdout_sha256": local["stdout_sha256"],
        },
    }


def incident_evidence():
    runs = []
    for _ in range(2):
        run = run_command([sys.executable, "checks.py", "incident"], ROOT)
        lines = run["stderr"].strip().splitlines()
        last = lines[-1] if lines else ""
        runs.append({
            "exit_code": run["exit_code"],
            "marker_present": MARKER in run["stderr"],
            "stdout_sha256": run["stdout_sha256"],
            # hash only the final line: tracebacks embed absolute paths
            "stderr_last_line_sha256": hashlib.sha256(last.encode()).hexdigest(),
            "stderr_last_line": last,
            "stdout_json_valid": _json_valid(run["stdout"]),
        })
    return runs


def _json_valid(text):
    try:
        json.loads(text)
        return True
    except (ValueError, TypeError):
        return False


def log_evidence():
    """Compare the committed logs against the independent model."""
    logs = {}
    for name in ("normal", "local", "incident-1", "incident-2"):
        path = ROOT / "logs" / ("%s.log" % name)
        text = path.read_text()
        lines = text.strip().splitlines()
        entry = {"sha256": hashlib.sha256(text.encode()).hexdigest(),
                 "line_count": len(lines)}
        if name in ("incident-1", "incident-2"):
            entry["payload"] = json.loads(lines[0])
            entry["last_line"] = lines[-1]
        else:
            entry["payload"] = json.loads(lines[0])
        logs[name] = entry
    return logs


def pipeline_results():
    """Observed entry under test; used only as a model cross-check."""
    out = {}
    for name in PLANS:
        run = run_command([sys.executable, "pipeline.py", "data/incident.json", name], ROOT)
        payload = json.loads(run["stdout"])
        out[name] = {
            "exit_code": run["exit_code"],
            "result": payload["result"],
            "block_bases": [b["basis_windows"] for b in payload["trace"]["blocks"]],
            "block_trace_loadings": [b["trace_loadings"] for b in payload["trace"]["blocks"]],
            "merged_information": payload["trace"]["merged"]["information"],
            "merged_rhs": payload["trace"]["merged"]["rhs"],
            "merged_rank": payload["trace"]["merged"]["information_rank"],
        }
    return out


def anchor_lines(canon, plans, breaks, contrasts):
    def d(x):
        return repr(float(x))

    def s(x):
        """Anchor text for a Fraction or an already-encoded [numerator, denominator]."""
        if isinstance(x, list):
            return "%d/%d" % (x[0], x[1])
        return "%d/%d" % (x.numerator, x.denominator)

    rows = [
        ("canonical.intercept", s(canon["beta"][0][0])),
        ("canonical.slope", s(canon["beta"][1][0])),
        ("canonical.slope_variance", s(canon["result"]["slope_variance"])),
        ("canonical.standard_error", d(canon["result"]["standard_error"])),
        ("canonical.ci95_low", d(canon["result"]["ci95"][0])),
        ("canonical.ci95_high", d(canon["result"]["ci95"][1])),
        ("canonical.information_rank", str(canon["result"]["information_rank"])),
        ("canonical.full_covariance_rank", str(canon["full_covariance_rank"])),
        ("canonical.intercept_variance", s(canon["covariance"][0][0])),
        ("canonical.basis", " ".join(canon["basis"])),
    ]
    for name in PLANS:
        p = plans[name]
        rows.extend([
            ("plans.%s.slope" % name, s(p["slope"])),
            ("plans.%s.reported_slope_variance" % name, s(p["reported_slope_variance"])),
            ("plans.%s.actual_slope_variance" % name, s(p["actual_slope_variance"])),
            ("plans.%s.variance_ratio" % name, s(p["variance_ratio"])),
            ("plans.%s.information_rank" % name, str(p["information_rank"])),
            ("plans.%s.union_window_rank" % name, str(p["union_window_rank"])),
            ("plans.%s.reported_intercept_variance" % name,
             s(p["merged_information_inverse"][0][0])),
            ("plans.%s.actual_intercept_variance" % name, s(p["actual_intercept_variance"])),
            ("plans.%s.ci95_low" % name, d(p["result"]["ci95"][0])),
            ("plans.%s.ci95_high" % name, d(p["result"]["ci95"][1])),
            ("plans.%s.actual_ci95_low" % name, d(p["actual_ci95"][0])),
            ("plans.%s.actual_ci95_high" % name, d(p["actual_ci95"][1])),
        ])
    dup = contrasts["duplicate_shared_source"]
    rows.extend([
        ("contrast.duplicate.canonical_slope", s(dup["canonical"]["slope"])),
        ("contrast.duplicate.canonical_slope_variance",
         s(dup["canonical"]["result"]["slope_variance"])),
        ("contrast.duplicate.canonical_rank", str(dup["canonical"]["result"]["information_rank"])),
        ("contrast.duplicate.all_slope", s(dup["all"]["slope"])),
        ("contrast.duplicate.all_reported_variance", s(dup["all"]["reported_slope_variance"])),
        ("contrast.duplicate.all_actual_variance", s(dup["all"]["actual_slope_variance"])),
        ("contrast.duplicate.all_rank", str(dup["all"]["information_rank"])),
        ("contrast.duplicate.singles_slope", s(dup["singles"]["slope"])),
        ("contrast.duplicate.singles_reported_variance",
         s(dup["singles"]["reported_slope_variance"])),
        ("contrast.duplicate.singles_actual_variance",
         s(dup["singles"]["actual_slope_variance"])),
        ("contrast.duplicate.singles_variance_ratio", s(dup["singles"]["variance_ratio"])),
        ("contrast.duplicate.singles_rank", str(dup["singles"]["information_rank"])),
    ])
    ind = contrasts["independent_sources"]
    rows.extend([
        ("contrast.independent.canonical_slope", s(ind["canonical"]["slope"])),
        ("contrast.independent.canonical_slope_variance",
         s(ind["canonical"]["result"]["slope_variance"])),
        ("contrast.independent.canonical_rank", str(ind["canonical"]["result"]["information_rank"])),
        ("contrast.independent.all_plans_equal_canonical", str(ind["all_plans_equal_canonical"])),
        ("contrast.independent.all_plans_ratio_one", str(ind["all_plans_ratio_one"])),
    ])
    for name in PLANS:
        d = breaks[name]
        rows.extend([
            ("decomposition.%s.rank_deficit" % name, str(d["rank_deficit"])),
            ("decomposition.%s.information_difference_rank" % name,
             str(d["information_difference_rank"])),
            ("decomposition.%s.information_difference_trace" % name,
             s(d["information_difference_trace"])),
            ("decomposition.%s.information_difference_determinant" % name,
             s(d["information_difference_determinant"])),
            ("decomposition.%s.information_difference_psd" % name,
             str(d["information_difference_psd"])),
            ("decomposition.%s.shift_matches" % name, str(d["shift_matches"])),
        ])
    rows.append(("decomposition.cut.phantom_is_source_c_squared",
                 str(breaks["cut"]["phantom_is_source_c_squared"])))
    return rows


def normal_evidence():
    """Independent recalculation of the healthy package (no shared sources)."""
    pkg = read_json(ROOT / "data" / "normal.json")
    batch = Batch(pkg, "normal")
    canon = batch.canonical()
    out = {
        "canonical": canonical_evidence(canon),
        "plans": {},
        "all_plans_invariant": True,
    }
    for name in sorted(batch.plans):
        ex = batch.existing(batch.plans[name])
        same = (ex["beta"] == canon["beta"]
                and ex["reported"]["slope_variance"] == canon["result"]["slope_variance"]
                and ex["actual_slope_variance"] == ex["reported"]["slope_variance"])
        out["all_plans_invariant"] = out["all_plans_invariant"] and same
        out["plans"][name] = {
            "slope": frac(ex["beta"][1][0]),
            "intercept": frac(ex["beta"][0][0]),
            "reported_slope_variance": frac(ex["reported"]["slope_variance"]),
            "actual_slope_variance": frac(ex["actual_slope_variance"]),
            "information_rank": ex["reported"]["information_rank"],
            "matches_canonical": same,
        }
    return out


def build_evidence():
    incident = read_json(ROOT / "data" / "incident.json")
    batch = Batch(incident, "incident")
    canon = batch.canonical()
    plans = plan_evidence(batch, canon)

    dup_batch = Batch(duplicate_package(incident), "duplicate")
    dup_canon = dup_batch.canonical()
    dup_plans = [["w0", "w1", "w2", "w3", "w4", "w5", "w6", "w7", "w8"]]
    dup_singles = [[w] for w in dup_batch.order]
    dup_all = dup_batch.existing(dup_plans)
    dup_single = dup_batch.existing(dup_singles)

    ind_batch = Batch(independent_package(incident), "independent")
    ind_canon = ind_batch.canonical()
    ind_plans = {}
    for name in PLANS:
        ex = ind_batch.existing(ind_batch.plans[name])
        ind_plans[name] = {
            "slope": frac(ex["beta"][1][0]),
            "reported_slope_variance": frac(ex["reported"]["slope_variance"]),
            "actual_slope_variance": frac(ex["actual_slope_variance"]),
            "information_rank": ex["reported"]["information_rank"],
            "equals_canonical_point": ex["beta"] == ind_canon["beta"],
            "equals_canonical_variance": ex["reported"]["slope_variance"] == ind_canon["result"]["slope_variance"],
            "ratio_is_one": ex["actual_slope_variance"] == ex["reported"]["slope_variance"],
        }

    contrasts = {
        "partition_only": {
            "description": "original package, only the block partition changes",
            "plans": {name: {
                "slope": plans[name]["slope"],
                "reported_slope_variance": plans[name]["reported_slope_variance"],
                "actual_slope_variance": plans[name]["actual_slope_variance"],
                "information_rank": plans[name]["information_rank"],
            } for name in PLANS},
        },
        "duplicate_shared_source": {
            "description": "append w8 with the same weights as w3 (b:c = 1/2:1/2)",
            "added_window": "w8",
            "canonical": canonical_evidence(dup_canon),
            "all": {
                "basis": [b["basis"] for b in dup_all["blocks"]],
                "slope": frac(dup_all["beta"][1][0]),
                "reported_slope_variance": frac(dup_all["reported"]["slope_variance"]),
                "actual_slope_variance": frac(dup_all["actual_slope_variance"]),
                "reported_information_equals_marginal_sum":
                    dup_all["reported_information_equals_marginal_sum"],
                "union_window_rank": dup_all["union_window_rank"],
                "union_nullspace": dup_all["union_nullspace"],
                "information_rank": dup_all["reported"]["information_rank"],
            },
            "singles": {
                "basis": [b["basis"] for b in dup_single["blocks"]],
                "slope": frac(dup_single["beta"][1][0]),
                "reported_slope_variance": frac(dup_single["reported"]["slope_variance"]),
                "actual_slope_variance": frac(dup_single["actual_slope_variance"]),
                "variance_ratio": frac(dup_single["variance_ratio"]),
                "reported_information_equals_marginal_sum":
                    dup_single["reported_information_equals_marginal_sum"],
                "union_window_rank": dup_single["union_window_rank"],
                "union_nullspace": dup_single["union_nullspace"],
                "information_rank": dup_single["reported"]["information_rank"],
            },
        },
        "independent_sources": {
            "description": "every window's referenced source replaced by a fresh id; "
                           "same time, value, variance and weights",
            "new_source_count": len(ind_batch.ids),
            "canonical": canonical_evidence(ind_canon),
            "plans": ind_plans,
            "all_plans_equal_canonical": all(
                v["equals_canonical_point"] and v["equals_canonical_variance"]
                for v in ind_plans.values()),
            "all_plans_ratio_one": all(v["ratio_is_one"] for v in ind_plans.values()),
        },
    }

    evidence = {
        "version": 1,
        "tolerance": TOL,
        "package_sha256": {
            "data/incident.json": hashlib.sha256(
                (ROOT / "data" / "incident.json").read_bytes()).hexdigest(),
            "data/normal.json": hashlib.sha256(
                (ROOT / "data" / "normal.json").read_bytes()).hexdigest(),
        },
        "full_batch": {
            "source_ids": batch.ids,
            "window_ids": batch.order,
            "window_covariance": encode(batch.full_cov()),
            "covariance_rank": rank(batch.full_cov()),
            "window_count": len(batch.windows),
        },
        "canonical": canonical_evidence(canon),
        "plans": plans,
        "decomposition": decomposition_evidence(batch, canon),
        "contrasts": contrasts,
        "logs": log_evidence(),
        "normal_model": normal_evidence(),
        "regressions": regression_evidence(),
        "incident_runs": incident_evidence(),
        "pipeline_cross_check": pipeline_results(),
    }
    evidence["anchors"] = [{"key": k, "value": v}
                           for k, v in anchor_lines(canon, plans,
                                                    evidence["decomposition"], contrasts)]
    return evidence


# --------------------------------------------------------------------------
# verify
# --------------------------------------------------------------------------
class VerifyError(AssertionError):
    pass


def need(condition, message):
    if not condition:
        raise VerifyError(message)


def close(a, b, tol=TOL):
    return abs(float(a) - float(b)) <= tol


def check_same(actual, stored, path):
    """Exact structural comparison of recomputed evidence with evidence.json."""
    need(type(actual) is type(stored) or (
        isinstance(actual, (int, float)) and isinstance(stored, (int, float))),
        "evidence type mismatch at %s: %r vs %r" % (path, actual, stored))
    if isinstance(stored, dict):
        need(set(actual) == set(stored),
             "evidence keys mismatch at %s: %r" % (path, sorted(set(actual) ^ set(stored))))
        for key in stored:
            check_same(actual[key], stored[key], "%s.%s" % (path, key))
    elif isinstance(stored, list):
        need(len(actual) == len(stored),
             "evidence length mismatch at %s: %d vs %d" % (path, len(actual), len(stored)))
        for i, (a, b) in enumerate(zip(actual, stored)):
            check_same(a, b, "%s[%d]" % (path, i))
    elif isinstance(stored, float):
        need(close(actual, stored), "evidence float mismatch at %s: %r vs %r"
             % (path, actual, stored))
    else:
        need(actual == stored, "evidence value mismatch at %s: %r vs %r"
             % (path, actual, stored))


def verify():
    evidence_path = ROOT / "evidence.json"
    report_path = ROOT / "REPORT.md"
    need(evidence_path.exists(), "evidence.json missing")
    need(report_path.exists(), "REPORT.md missing")
    stored = json.loads(evidence_path.read_text())
    fresh = build_evidence()

    # 1. recomputed evidence must match the committed evidence
    check_same(fresh, stored, "evidence")

    # 2. numeric anchors must appear in REPORT.md
    report = report_path.read_text()
    for anchor in fresh["anchors"]:
        line = "%s = %s" % (anchor["key"], anchor["value"])
        need(line in report, "REPORT.md missing anchor: %s" % line)

    # 3. original normal / local regressions
    norm = fresh["regressions"]["normal"]
    need(norm["exit_code"] == 0, "checks.py normal failed")
    for name, result in norm["result"].items():
        need(result["intercept"] == [1, 1], "normal intercept drift (%s)" % name)
        need(result["slope"] == [2, 1], "normal slope drift (%s)" % name)
        need(result["slope_variance"] == [1, 5], "normal variance drift (%s)" % name)
        need(result["information_rank"] == 4, "normal rank drift (%s)" % name)
    need(norm["result"]["all"] == norm["result"]["cut"] == norm["result"]["single"],
         "normal plan invariance broken")
    nm = fresh["normal_model"]
    need(nm["all_plans_invariant"], "independent model disagrees across normal plans")
    for name, model in nm["plans"].items():
        need(model["intercept"] == [1, 1], "independent normal intercept (%s)" % name)
        need(model["slope"] == [2, 1], "independent normal slope (%s)" % name)
        need(model["reported_slope_variance"] == [1, 5],
             "independent normal variance (%s)" % name)
        need(model["actual_slope_variance"] == [1, 5],
             "independent normal actual variance (%s)" % name)
        need(model["information_rank"] == 4, "independent normal rank (%s)" % name)
        need(model["matches_canonical"], "independent normal canonical mismatch (%s)" % name)
    loc = fresh["regressions"]["local"]
    need(loc["exit_code"] == 0, "checks.py local failed")
    need(loc["result"].get("local_covariance_checks") == "passed",
         "local covariance checks not passed")

    # 4. the original incident entry must fail twice in a row, same marker
    runs = fresh["incident_runs"]
    need(len(runs) == 2, "expected two incident runs")
    for i, run in enumerate(runs):
        need(run["exit_code"] != 0, "incident run %d unexpectedly succeeded" % (i + 1))
        need(run["marker_present"], "incident run %d missing %s" % (i + 1, MARKER))
        need(run["stdout_json_valid"], "incident run %d stdout is not JSON" % (i + 1))
    need(runs[0]["stdout_sha256"] == runs[1]["stdout_sha256"],
         "incident stdout differs between consecutive runs")
    need(runs[0]["stderr_last_line_sha256"] == runs[1]["stderr_last_line_sha256"],
         "incident stderr differs between consecutive runs")

    # 5. the alive logs must agree with the recomputed model
    logs = fresh["logs"]
    for name, result in logs["normal"]["payload"]["normal"].items():
        need(result == fresh["regressions"]["normal"]["result"][name],
             "logs/normal.log disagrees with replayed normal plan %s" % name)
    need(logs["local"]["payload"] == fresh["regressions"]["local"]["result"],
         "logs/local.log disagrees with replayed local checks")
    for name in ("incident-1", "incident-2"):
        entry = logs[name]
        need(entry["last_line"] == "AssertionError: " + MARKER,
             "logs/%s.log does not end with the incident marker" % name)
        for plan, observed in entry["payload"]["incident"].items():
            model = fresh["plans"][plan]["result"]
            for key in ("intercept", "slope", "slope_variance",
                        "information_rank", "window_count"):
                need(observed[key] == model[key],
                     "logs/%s.log %s %s disagrees with model" % (name, plan, key))
            need(close(observed["standard_error"], model["standard_error"]),
                 "logs/%s.log %s standard_error disagrees" % (name, plan))
            for a, b in zip(observed["ci95"], model["ci95"]):
                need(close(a, b), "logs/%s.log %s ci95 disagrees" % (name, plan))

    # 6. the independent model must reproduce what the original entry reports
    batch = Batch(read_json(ROOT / "data" / "incident.json"), "incident")
    for name, observed in fresh["pipeline_cross_check"].items():
        need(observed["exit_code"] == 0, "pipeline.py %s failed" % name)
        result = observed["result"]
        model = fresh["plans"][name]["result"]
        for key in ("information_rank", "window_count"):
            need(result[key] == model[key],
                 "pipeline %s %s mismatch: %r vs %r"
                 % (name, key, result[key], model[key]))
        for key in ("intercept", "slope", "slope_variance"):
            need(result[key] == model[key],
                 "pipeline %s %s mismatch: %r vs %r"
                 % (name, key, result[key], model[key]))
        need(close(result["standard_error"], model["standard_error"]),
             "pipeline %s standard_error mismatch" % name)
        for a, b in zip(result["ci95"], model["ci95"]):
            need(close(a, b), "pipeline %s ci95 mismatch" % name)

        # window lineage and merged score must equal the independent model
        model_blocks = batch.existing_blocks(batch.plans[name])
        need(observed["block_bases"] == [b["basis"] for b in model_blocks],
             "pipeline %s basis mismatch" % name)
        for got, block in zip(observed["block_trace_loadings"], model_blocks):
            want = {}
            for k, (c0, c1) in zip(batch.ids, zip(block["loadings"][0],
                                                  block["loadings"][1])):
                if c0 or c1:
                    want[k] = [frac(c0), frac(c1)]
            need({k: [[q[0], q[1]] for q in v] for k, v in got.items()} == want,
                 "pipeline %s trace_loadings mismatch" % name)
        need(observed["merged_information"] ==
             encode(fresh["plans"][name]["merged_information"]),
             "pipeline %s merged information mismatch" % name)
        need(observed["merged_rhs"] ==
             encode(fresh["plans"][name]["merged_rhs"]),
             "pipeline %s merged rhs mismatch" % name)
        need(observed["merged_rank"] == fresh["plans"][name]["information_rank"],
             "pipeline %s merged rank mismatch" % name)

    print(json.dumps({
        "verify": "ok",
        "evidence": "matches recomputation",
        "report_anchors": len(fresh["anchors"]),
        "normal_exit": norm["exit_code"],
        "local_exit": loc["exit_code"],
        "incident_exits": [r["exit_code"] for r in runs],
        "incident_marker": MARKER,
        "plans_cross_checked": sorted(fresh["plans"]),
    }, sort_keys=True))


# --------------------------------------------------------------------------
# demo
# --------------------------------------------------------------------------
def demo():
    incident = read_json(ROOT / "data" / "incident.json")
    batch = Batch(incident, "incident")
    canon = batch.canonical()
    cut = batch.existing(batch.plans["cut"])
    dup_batch = Batch(duplicate_package(incident), "duplicate")
    dup_singles = dup_batch.existing([[w] for w in dup_batch.order])

    payload = {
        "canonical": {
            "basis": canon["basis"],
            "intercept": frac(canon["beta"][0][0]),
            "slope": frac(canon["beta"][1][0]),
            "slope_variance": frac(canon["result"]["slope_variance"]),
            "standard_error": canon["result"]["standard_error"],
            "ci95": canon["result"]["ci95"],
            "information_rank": canon["result"]["information_rank"],
        },
        "existing_cut": {
            "block_bases": [b["basis"] for b in cut["blocks"]],
            "slope": frac(cut["beta"][1][0]),
            "reported_slope_variance": frac(cut["reported"]["slope_variance"]),
            "reported_ci95": cut["reported"]["ci95"],
            "actual_slope_variance": frac(cut["actual_slope_variance"]),
            "actual_ci95": cut["actual_ci95"],
            "information_rank": cut["reported"]["information_rank"],
        },
        "single_factor_difference": {
            "factor": "append w8 duplicating w3's sources and weights",
            "canonical_slope_variance": frac(
                dup_batch.canonical()["result"]["slope_variance"]),
            "all_one_block_reported_variance": frac(
                dup_batch.existing([[w for w in dup_batch.order]])["reported"]["slope_variance"]),
            "all_single_blocks_reported_variance": frac(
                dup_singles["reported"]["slope_variance"]),
            "all_single_blocks_actual_variance": frac(dup_singles["actual_slope_variance"]),
            "all_single_blocks_information_rank": dup_singles["reported"]["information_rank"],
        },
    }
    print(json.dumps(encode(payload), ensure_ascii=False, sort_keys=True))


# --------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["build", "verify", "demo"])
    args = parser.parse_args()
    if args.mode == "build":
        evidence = build_evidence()
        print(json.dumps(encode(evidence), ensure_ascii=False, sort_keys=True, indent=2))
    elif args.mode == "verify":
        verify()
    else:
        demo()


if __name__ == "__main__":
    main()

