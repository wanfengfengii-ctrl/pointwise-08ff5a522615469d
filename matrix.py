from fractions import Fraction as F

def rank(rows):
    a = [[F(x) for x in row] for row in rows]
    if not a:
        return 0
    r = 0
    for c in range(len(a[0])):
        p = next((i for i in range(r, len(a)) if a[i][c]), None)
        if p is None:
            continue
        a[r], a[p] = a[p], a[r]
        d = a[r][c]
        a[r] = [x / d for x in a[r]]
        for i in range(r + 1, len(a)):
            d = a[i][c]
            a[i] = [x - d*y for x, y in zip(a[i], a[r])]
        r += 1
        if r == len(a):
            break
    return r

def unpack(q):
    return F(q[0], q[1]) if isinstance(q, list) else F(q)

def encode(x):
    if isinstance(x, F):
        return [x.numerator, x.denominator]
    if isinstance(x, dict):
        return {k: encode(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [encode(v) for v in x]
    return x

def transpose(a):
    return [list(c) for c in zip(*a)]

def multiply(a, b):
    return [[sum((x*y for x, y in zip(row, col)), F(0)) for col in zip(*b)] for row in a]

def inverse(a):
    n = len(a)
    m = [[F(x) for x in row] + [F(i == j) for j in range(n)] for i, row in enumerate(a)]
    for c in range(n):
        p = next((r for r in range(c, n) if m[r][c]), None)
        if p is None:
            raise ValueError("singular matrix")
        m[c], m[p] = m[p], m[c]
        d = m[c][c]
        m[c] = [x/d for x in m[c]]
        for r in range(n):
            if r != c:
                d = m[r][c]
                m[r] = [x-d*y for x, y in zip(m[r], m[c])]
    return [r[n:] for r in m]
