"""Exact scalar and small-matrix checks for the two ordinary-group addenda."""

from fractions import Fraction
import json
from math import comb, gcd, lcm, prod

from verify import cyclotomic_floor
from sympy import primerange


def remaining_product(a, b, c):
    return prod(
        Fraction(c + i + j - 1, i + j - 1)
        for i in range(2, a + 1) for j in range(1, b)
    )


def addendum_evaluation(b, c):
    upper = b + c
    lower = comb(upper, b)
    gcd_of_entries = gcd(comb(upper, b - 1), lower, comb(upper, b + 1))
    cleared_ratio = lcm(
        (c + 1) // gcd(c + 1, b),
        (b + 1) // gcd(b + 1, c),
    )
    expected_first = 1
    expected_second = 1
    for p in primerange(2, upper + 2):
        power = p
        while power < upper + 2:
            e = cyclotomic_floor(2, b, c, power)
            if e == 2:
                expected_first *= p
            if e:
                expected_second *= p
            power *= p
    determinant = lower**2 - comb(upper, b - 1) * comb(upper, b + 1)
    assert lower == cleared_ratio * gcd_of_entries
    assert expected_first == gcd_of_entries
    assert expected_first * expected_second == determinant
    return [expected_first, expected_second]


def main():
    table = {}
    for b in (4, 5, 6):
        f3, f4 = remaining_product(3, b, b), remaining_product(4, b, b)
        l3, l4 = lcm(*range(1, b + 3)), lcm(*range(1, 2 * b))
        assert f3 > l3 and f4 > l4
        table[b] = [str(f3), l3, str(f4), l4]
    assert remaining_product(3, 7, 7) == Fraction(306735, 4) > 4**8
    assert remaining_product(3, 3, 6) == 90 > 60
    assert [cyclotomic_floor(3, 3, c, 7) for c in (3, 4, 5)] == [2, 3, 2]
    assert [addendum_evaluation(2, c)[0] for c in (2, 3, 4)] == [2, 5, 1]
    assert addendum_evaluation(2, 4) == [1, 105]

    witnesses = [(1, 7), (5, 7), (7, 15), (30, 30), (5, 22), (8, 17)]
    smith = {f"2,{b},{c}": addendum_evaluation(b, c) for b, c in witnesses}
    levels = {
        d: cyclotomic_floor(2, 30, 30, d)
        for d in (2, 4, 8, 16, 32)
        if cyclotomic_floor(2, 30, 30, d)
    }
    assert levels == {4: 2, 8: 2, 16: 2, 32: 2}
    print(json.dumps({"status": "PASS", "cyclicity_table": table,
                      "side_two_smith_examples": smith,
                      "four_active_dyadic_levels": levels}, indent=2))


if __name__ == "__main__":
    main()
