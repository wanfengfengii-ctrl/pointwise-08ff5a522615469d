import argparse
import json
from fractions import Fraction as F
from pathlib import Path
from matrix import encode, inverse, multiply
from pipeline import run

MARKER = "TREND_PARTITION_DRIFT"

def load(path):
    return json.loads(Path(path).read_text())

def normal():
    p = load("data/normal.json")
    results = {name: run(p, plan)["result"] for name, plan in p["plans"].items()}
    for r in results.values():
        assert r["slope"] == 2 and r["intercept"] == 1
        assert r["slope_variance"] == F(1, 5)
        assert r["information_rank"] == 4
    assert all(r == results["all"] for r in results.values())
    print(json.dumps(encode({"normal": results}), sort_keys=True))

def local():
    count = 0
    for path in ("data/normal.json", "data/incident.json"):
        p = load(path)
        for plan in p["plans"].values():
            out = run(p, plan)
            for block in out["trace"]["blocks"]:
                c = block["covariance"]
                assert all(c[i][i] > 0 for i in range(len(c)))
                assert all(c[i][j] == c[j][i] for i in range(len(c)) for j in range(len(c)))
                unit = multiply(c, inverse(c))
                assert all(unit[i][j] == int(i == j) for i in range(len(c)) for j in range(len(c)))
                a = block["information"]
                assert a[0][1] == a[1][0] and a[0][0] > 0
                assert a[0][0]*a[1][1] - a[0][1]*a[1][0] >= 0
                count += 1
    print(json.dumps({"checked_block_count": count, "local_covariance_checks": "passed"}))

def incident():
    p = load("data/incident.json")
    rows = {name: run(p, plan)["result"] for name, plan in p["plans"].items()}
    print(json.dumps(encode({"incident": rows}), sort_keys=True), flush=True)
    signature = lambda r: (r["slope_variance"], r["information_rank"])
    if any(signature(r) != signature(rows["all"]) for r in rows.values()):
        raise AssertionError(MARKER)
    print("partition invariant")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["normal", "local", "incident"])
    args = parser.parse_args()
    globals()[args.mode]()

if __name__ == "__main__":
    main()
