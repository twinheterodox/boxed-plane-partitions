"""Verify every Laurent determinantal divisor on fixed boxes. Requires SymPy."""

import json
from itertools import combinations
from math import comb
from pathlib import Path

from sympy import QQ, ZZ, Symbol, cyclotomic_poly
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import invariant_factors


BOXES = (
    (1, 1, 1),
    (1, 4, 5),
    (3, 1, 2),
    (3, 2, 1),
    (2, 2, 2),
    (3, 4, 3),
    (4, 3, 5),
    (4, 4, 4),
    (5, 3, 3),
)
SYMBOL = Symbol("q")
RING = QQ.poly_ring(SYMBOL)
Q = RING.gens[0]


def carlitz_matrix(a, b, c):
    elementary = [RING.one]
    for j in range(b + c):
        elementary.append(RING.zero)
        for k in range(j + 1, 0, -1):
            elementary[k] += Q**j * elementary[k - 1]

    def entry(i, j):
        k = b - i + j
        return elementary[k] if 0 <= k < len(elementary) else RING.zero

    rows = [[entry(i, j) for j in range(a)] for i in range(a)]
    return DomainMatrix(rows, (a, a), RING)


def volume_multiplicities(a, b, c):
    multiplicities = {}
    for d in range(2, a + b + c):
        count = 0
        for i in range(1, b + 1):
            for j in range(1, c + 1):
                count += int((a + i + j - 1) % d == 0)
                count -= int((i + j - 1) % d == 0)
        if count:
            assert count > 0
            multiplicities[d] = count
    return multiplicities


def laurent_monic(poly):
    order_at_zero = min(exponents[0] for exponents in poly)
    return RING.exquo(poly, Q**order_at_zero).monic()


def check_box(a, b, c):
    matrix = carlitz_matrix(a, b, c)
    multiplicities = volume_multiplicities(a, b, c)
    factors = {
        d: RING.from_sympy(cyclotomic_poly(d, SYMBOL))
        for d in multiplicities
    }
    count = 0
    for size in range(1, a + 1):
        predicted = RING.one
        for d, e in multiplicities.items():
            predicted *= factors[d] ** max(0, e - a + size)
        indices = list(combinations(range(a), size))
        divisor = RING.zero
        for rows in indices:
            for columns in indices:
                minor = matrix.extract(rows, columns).det()
                divisor = RING.gcd(divisor, minor)
                count += 1
        assert laurent_monic(divisor) == predicted, (a, b, c, size)

    def integer_entry(i, j):
        k = b - i + j
        return ZZ(comb(b + c, k)) if 0 <= k <= b + c else ZZ.zero

    entries = [[integer_entry(i, j) for j in range(a)] for i in range(a)]
    integer_factors = invariant_factors(DomainMatrix(entries, (a, a), ZZ))
    return {
        "box": [a, b, c],
        "minors": count,
        "macmahon_cyclotomic_multiplicities": multiplicities,
        "ordinary_smith_factors": [abs(int(f)) for f in integer_factors],
    }


if __name__ == "__main__":
    boxes = [check_box(*box) for box in BOXES]
    result = {
        "method": "gcd of all minors of each size in QQ[q], discarding units of the Laurent ring",
        "boxes": boxes,
        "minors_total": sum(item["minors"] for item in boxes),
        "status": "PASS",
    }
    output = Path(__file__).with_name("smith-checks.json")
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"PASS on {len(boxes)} boxes and {result['minors_total']} polynomial minors")
