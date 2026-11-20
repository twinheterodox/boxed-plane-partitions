"""Test the interval-ideal lemma directly by polynomial gcds in finite characteristic."""

import json
from pathlib import Path

from sympy import GF, QQ, Symbol


Q_SYMBOL = Symbol("q")
PRIMES = (2, 3, 5, 7, 11)
MAX_ROW = 12


def gaussian_rows(ring):
    q = ring.gens[0]
    row = [ring.one]
    yield row
    for n in range(1, MAX_ROW + 1):
        row = [
            (row[k - 1] if k else ring.zero)
            + (q**k * row[k] if k < n else ring.zero)
            for k in range(n + 1)
        ]
        yield row


def interval_gcds(row, ring):
    for left in range(len(row)):
        divisor = ring.zero
        for right in range(left, len(row)):
            divisor = ring.gcd(divisor, row[right])
            yield left, right, divisor.monic()


def verify():
    rational_ring = QQ.poly_ring(Q_SYMBOL)
    exact_rows = list(gaussian_rows(rational_ring))
    counts = {}
    for p in PRIMES:
        finite_ring = GF(p).poly_ring(Q_SYMBOL)
        checks = 0
        for n, row in enumerate(gaussian_rows(finite_ring)):
            expected = interval_gcds(exact_rows[n], rational_ring)
            actual = interval_gcds(row, finite_ring)
            for (u, v, gcd_q), (x, y, gcd_p) in zip(expected, actual, strict=True):
                assert (u, v) == (x, y)
                assert gcd_p == finite_ring.from_sympy(gcd_q.as_expr()), (p, n, u, v)
                checks += 1
        counts[p] = checks
    return {
        "method": "compute the whole polynomial gcd before and after reduction mod p",
        "all_nonempty_intervals_in_rows": [0, MAX_ROW],
        "intervals_per_prime": counts,
        "total": sum(counts.values()),
        "status": "PASS",
    }


if __name__ == "__main__":
    result = verify()
    output = Path(__file__).with_name("interval-checks.json")
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"PASS on {result['total']} interval/prime pairs")
