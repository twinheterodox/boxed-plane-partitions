"""Compare the arithmetic formulas with ordinary Smith forms and Lucas ranks."""

import json
from collections import Counter
from functools import cache
from itertools import product
from math import comb, prod

import sympy
from sympy.matrices.normalforms import smith_normal_form


def cyclotomic_floor(a, b, c, d):
    return sum(
        (b + c + i) // d + i // d - (b + i) // d - (c + i) // d
        for i in range(a)
    )


def cyclotomic_residue(a, b, c, d):
    x, y, z = (side % d for side in (a, b, c))
    s = x + y + z
    return max(0, min(x, y, z, d - x, d - y, d - z, s - d, 2 * d - s))


def lucas_rank(exponents, height, p):
    @cache
    def rank(columns, rows):
        if not columns or not rows:
            return 0
        if rows == 1:
            return 1
        full, partial = divmod(rows, p)
        slices = {}
        for e in columns:
            digit = e % p
            slices.setdefault(digit, []).append(e // p)
        result = 0
        active = 0
        for values in slices.values():
            lower = tuple(values)
            first = rank(lower, full)
            result += first
            if partial and rank(lower, full + 1) > first:
                active += 1
        return result + min(partial, active)

    return rank(tuple(exponents), height)


def binomial_or_zero(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def valuation(m, p):
    count = 0
    while m % p == 0:
        m //= p
        count += 1
    return count


def check_box(a, b, c, stats):
    total = a + b + c
    matrix = sympy.Matrix(
        a, a, lambda i, j: binomial_or_zero(b + c, b + i - j)
    )
    smith = smith_normal_form(matrix, domain=sympy.ZZ)
    diagonal = [abs(int(smith[i, i])) for i in range(a)]
    assert all(diagonal)
    assert abs(int(matrix.det())) == prod(diagonal)

    coefficients = {}
    for d in range(2, total + 2):
        e = cyclotomic_floor(a, b, c, d)
        assert cyclotomic_residue(a, b, c, d) == e
        coefficients[d] = e
        stats["floor_vs_residue"] += 1

    columns = list(range(b)) + list(range(a + b, total))
    report = {}
    for p in sympy.primerange(2, total + 2):
        p = int(p)
        d = p
        levels = []
        while d < total:
            if coefficients[d]:
                levels.append((d, coefficients[d]))
            d *= p
        counts = [e for _, e in levels]
        assert counts == sorted(counts)
        primary = [valuation(x, p) for x in diagonal if x % p == 0]
        nullity = max(counts, default=0)
        assert nullity == len(primary)
        assert nullity == len(columns) - lucas_rank(columns, b + c, p)
        assert sum(counts) == sum(primary)
        assert max(primary, default=0) <= len(counts)
        stats["rank_length_exponent_and_pascal_rank"] += 1

        if len(counts) <= 2:
            u, v = [0] * (2 - len(counts)) + counts
            assert primary == [1] * (v - u) + [2] * u
            stats["primary_groups_at_most_two_levels"] += 1
            if len(counts) == 2:
                stats["exactly_two_nonempty_levels"] += 1
                if levels[1][0] != p * levels[0][0]:
                    stats["nonconsecutive_two_levels"] += 1
        else:
            stats["three_or_more_nonempty_levels"] += 1
        if primary:
            report[p] = dict(levels=levels, smith_valuations=primary)
    stats["boxes"] += 1
    return dict(sides=(a, b, c), smith=diagonal, primary=report)


def main():
    boxes = set(product(range(1, 8), repeat=3))
    boxes.update([(3, 5, 5), (5, 8, 12), (8, 7, 9), (4, 16, 23)])
    examples = []
    stats = Counter()
    for sides in sorted(boxes):
        checked = check_box(*sides, stats)
        if sides in [(2, 2, 2), (3, 3, 3), (3, 5, 5)]:
            examples.append(checked)
    print(json.dumps(dict(status="PASS", sympy=sympy.__version__,
                          checks=dict(stats), examples=examples), indent=2))


if __name__ == "__main__":
    main()
