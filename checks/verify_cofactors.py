"""Exact polynomial checks of the boxed-plane-partition cofactor proof."""

from functools import cache
from math import comb, prod

from sympy.polys.domains import ZZ
from sympy.polys.rings import ring


R, q = ring("q", ZZ)


@cache
def cyclotomic(d):
    proper = (cyclotomic(j) for j in range(1, d) if d % j == 0)
    return (q**d - 1).exquo(prod(proper, start=R.one))


def elementary_functions(m):
    values = [R.one] + [R.zero for _ in range(m)]
    for j in range(m):
        for k in range(j + 1, 0, -1):
            values[k] = values[k] + q**j * values[k - 1]
    return values


@cache
def straight_schur(shape, m):
    padded = shape + (0,) * (m - len(shape))
    assert len(padded) == m
    numerator = q**sum(i * part for i, part in enumerate(padded))
    denominator = R.one
    for i in range(m):
        for j in range(i + 1, m):
            numerator *= 1 - q**(padded[i] - padded[j] + j - i)
            denominator *= 1 - q**(j - i)
    return numerator.exquo(denominator)


def strip_power_of_q(poly):
    power = min(monomial[0] for monomial in poly)
    return poly.exquo(q**power)


def check_box(a, b, c):
    m = b + c
    values = elementary_functions(m)
    matrix = [
        [values[h] if 0 <= h <= m else R.zero for h in range(b-i, b-i+a)]
        for i in range(a)
    ]
    for i in range(a):
        for j in range(a):
            h = b - i + j
            expected = comb(m, h) if 0 <= h <= m else 0
            assert matrix[i][j].evaluate(q, 1) == expected

    @cache
    def determinant(rows, columns):
        if not rows:
            return R.one
        total = R.zero
        for i, col in enumerate(columns):
            remaining = columns[:i] + columns[i+1:]
            total += (-1)**i * matrix[rows[0]][col] * determinant(rows[1:], remaining)
        return total

    indices = tuple(range(a))
    arms = [i+j-1 for i in range(1, b+1) for j in range(1, c+1)]
    multiplicities = {}
    for d in range(2, a+b+c):
        order = sum(int((a+k) % d == 0) - int(k % d == 0) for k in arms)
        assert 0 <= order <= a
        if order:
            multiplicities[d] = order
    macmahon = prod((cyclotomic(d)**e for d, e in multiplicities.items()), start=R.one)
    divisor = prod((cyclotomic(d)**(e-1) for d, e in multiplicities.items()), start=R.one)
    radical = macmahon.exquo(divisor)
    assert determinant(indices, indices) == q**(a*comb(b, 2)) * macmahon

    gcd = R.zero
    for r in indices:
        rows = indices[:r] + indices[r+1:]
        for s in indices:
            columns = indices[:s] + indices[s+1:]
            minor = determinant(rows, columns)
            # r, s are zero-based here, so the two partitions use k=r, l=s.
            pieri = R.zero
            for z in range(r, a):
                y = a - 1 + r - s - z
                if 0 <= y <= r:
                    shape = (a-1,)*(b-1) + (z, y)
                    pieri += straight_schur(shape, m)
            assert minor == pieri, (a, b, c, r, s)
            minor.exquo(divisor)
            assert (minor.evaluate(q, 1) * radical.evaluate(q, 1)) % macmahon.evaluate(q, 1) == 0
            gcd = gcd.gcd(minor)
    assert strip_power_of_q(gcd) == divisor
    return a*a, max(multiplicities.values(), default=0)


if __name__ == "__main__":
    boxes = (
        (1, 1, 1), (1, 3, 4), (2, 2, 2), (3, 1, 4),
        (3, 4, 1), (3, 2, 4), (4, 3, 2), (4, 4, 4),
    )
    total = 0
    for box in boxes:
        count, largest = check_box(*box)
        print(f"{box}  {count} minors  largest cyclotomic multiplicity {largest}")
        total += count
    print(f"PASS  {len(boxes)} boxes, {total} unsigned cofactors in Z[q]")
