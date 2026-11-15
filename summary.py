from fractions import Fraction as F
from matrix import rank, transpose, multiply, inverse

def summarize(source, windows):
    ids = sorted(source)
    rows, selected = [], []
    for w in windows:
        row = [w["weights"].get(k, F(0)) for k in ids]
        if rank(rows + [row]) > len(rows):
            rows.append(row)
            selected.append(w)
    covariance = [[sum((u*v*source[k]["variance"] for k, u, v in zip(ids, a, b)), F(0))
                   for b in rows] for a in rows]
    precision = inverse(covariance)
    design = [[F(1), w["time"]] for w in selected]
    observed = [[w["value"]] for w in selected]
    score = multiply(transpose(design), precision)
    information = multiply(score, design)
    rhs = multiply(score, observed)
    loadings = multiply(score, rows)
    return {
        "windows": [w["id"] for w in windows],
        "basis_windows": [w["id"] for w in selected],
        "covariance": covariance, "information": information, "rhs": rhs,
        "information_rank": len(rows),
        "trace_loadings": {k: [loadings[0][i], loadings[1][i]] for i, k in enumerate(ids)
                           if loadings[0][i] or loadings[1][i]}
    }

def merge(blocks):
    return {
        "information": [[sum(b["information"][i][j] for b in blocks) for j in range(2)] for i in range(2)],
        "rhs": [[sum(b["rhs"][i][0] for b in blocks)] for i in range(2)],
        "information_rank": sum(b["information_rank"] for b in blocks),
        "window_count": sum(len(b["windows"]) for b in blocks)
    }
