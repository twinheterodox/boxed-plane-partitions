"""Exact spot checks of the 19 November boxed plane partitions PDF. Python and SymPy.

Run this file directly. It has no manuscript or companion-checker dependencies.
"""
from collections import Counter
from fractions import Fraction
from functools import cache
from itertools import combinations, combinations_with_replacement, product
from math import comb, gcd, lcm, prod

from sympy import GF, ZZ, Matrix, Poly, Symbol, cyclotomic_poly, primerange
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.rings import ring


Zq, q = ring("q", ZZ)
symbol = Symbol("q")


@cache
def phi(d):
    return Zq.from_dict(Poly(cyclotomic_poly(d, symbol), symbol).as_dict())


def macmahon(a, b, c):
    return prod(Fraction(a+i+j-1, i+j-1)
                for i in range(1, b+1) for j in range(1, c+1))


def multiplicities(a, b, c):
    """Count divisible numerator and denominator indices in MacMahon's product."""
    indices = Counter()
    for i in range(1, b+1):
        for j in range(1, c+1):
            indices[a+i+j-1] += 1
            indices[i+j-1] -= 1
    return {d: sum(indices[k] for k in range(d, a+b+c, d))
            for d in range(2, a+b+c)}


def tent(a, b, c, d):
    r = [a % d, b % d, c % d]
    return max(0, min(*r, *(d-v for v in r), sum(r)-d, 2*d-sum(r)))


def elementary(n):
    coefficients = [Zq.one]
    for j in range(n):
        old = coefficients
        coefficients = old + [Zq.zero]
        for k, p in enumerate(old, start=1):
            coefficients[k] += q**j * p
    return coefficients


def strip_unit(p):
    primitive = p.primitive()[1]
    low = min(i for (i,) in primitive)
    result = primitive.exquo(q**low)
    return result if result.LC > 0 else -result


def ordinary(a, b, c):
    def entry(i, j):
        k = b-i+j
        return comb(b+c, k) if 0 <= k <= b+c else 0
    return Matrix(a, a, entry)


def integer_factors(matrix):
    smith = smith_normal_form(matrix, domain=ZZ)
    return tuple(abs(int(smith[i, i])) for i in range(smith.rows))


def polynomial_minors(a, b, c):
    es = elementary(b+c)
    entries = [[es[b-i+j] if 0 <= b-i+j <= b+c else Zq.zero
                for j in range(a)] for i in range(a)]
    actual_integers = integer_factors(ordinary(a,b,c))
    orders = multiplicities(a,b,c)

    @cache
    def det(rows, cols):
        if not rows:
            return Zq.one
        return sum(((-1)**j * entries[rows[0]][col]
                    * det(rows[1:], cols[:j]+cols[j+1:])
                    for j, col in enumerate(cols)), Zq.zero)

    total = 0
    for size in range(1, a+1):
        subsets = list(combinations(range(a), size))
        poly_gcd, int_gcd = Zq.zero, 0
        for rows, cols in product(subsets, repeat=2):
            value = det(rows, cols)
            poly_gcd = poly_gcd.gcd(value)
            int_gcd = gcd(int_gcd, int(value.evaluate(q, 1)))
            total += 1
        expected = prod((phi(d)**(e-a+size) for d,e in orders.items()
                         if e > a-size), start=Zq.one)
        assert strip_unit(poly_gcd) == expected, (a,b,c,size)
        assert int_gcd == prod(actual_integers[:size])
        if size == a:
            assert value == q**(a*comb(b,2)) * expected
    return total, actual_integers


def primary_checks(a,b,c):
    integers = integer_factors(ordinary(a,b,c))
    order = prod(integers)
    assert order == macmahon(a,b,c)
    assert (integers[-1] == order) == (min(a,b,c)==1 or sorted((a,b,c))==[2,2,4])
    exps = multiplicities(a,b,c)
    assert exps == {d: tent(a,b,c,d) for d in exps}
    assert max(exps.values()) == min(a,b,c)
    primes, classified = 0, 0
    predicted_first = 1
    for p in primerange(2,a+b+c):
        multiplicities_p = []
        power = p
        while power < a+b+c:
            if exps[power]:
                multiplicities_p.append(exps[power])
            if exps[power] == a:
                predicted_first *= p
            power *= p
        vals = []
        for number in integers:
            v = 0
            while number % p == 0:
                number //= p
                v += 1
            if v:
                vals.append(v)
        assert len(vals) == max(multiplicities_p,default=0)
        assert sum(vals) == sum(multiplicities_p)
        assert max(vals,default=0) <= len(multiplicities_p)
        if len(multiplicities_p) <= 2:
            u,v = sorted((multiplicities_p + [0,0])[:2])
            assert vals == [1]*(v-u) + [2]*u
            classified += 1
        primes += 1
    assert integers[0] == predicted_first
    return primes, classified


def interval_checks(limit=12):
    """Form Gaussian polynomials by the q-Pascal recurrence, then reduce coefficients."""
    fields = [ring("q", GF(p))[0] for p in (2,3,5,7,11)]
    row = [Zq.one]
    tested = 0
    for n in range(limit+1):
        for start in range(n+1):
            common = Zq.zero
            mod_gcds = [field.zero for field in fields]
            for entry in row[start:]:
                common = common.gcd(entry)
                for i, field in enumerate(fields):
                    mod_gcds[i] = mod_gcds[i].gcd(field.from_dict(dict(entry)))
                    assert mod_gcds[i].monic() == field.from_dict(dict(common)).monic()
                    tested += 1
        row = [Zq.one] + [row[k-1] + q**k*row[k] for k in range(1,n+1)] + [Zq.one]
    return tested


def check_examples():
    volumes = Counter(sum(h) for h in product(range(3), repeat=4)
                      if h[0]>=h[1]>=h[3] and h[0]>=h[2]>=h[3])
    assert [volumes[i] for i in range(9)] == [1,1,3,3,4,3,3,1,1]
    assert macmahon(1,3,3) == macmahon(2,2,2) == 20
    assert macmahon(2,2,4) == 105
    c = ordinary(2,2,4)
    assert all(v % 105 == 0 for v in Matrix([[-9,5]])*c)
    assert -9+2*5 == 1
    assert integer_factors(ordinary(3,3,3)) == (1,7,140)
    assert integer_factors(ordinary(9,9,9))[-1] % 4 == 0
    for h in range(2,9):
        k=2**h+1
        powers={2**j:tent(k,k,k,2**j) for j in range(1,(3*k-1).bit_length()+1)}
        assert {j:e for j,e in powers.items() if e} == {2:1,2**(h+1):2**h-3}
    for a,b,c in [(3,3,3),(3,3,4),(3,3,5)]:
        print("small cube family", (a,b,c), macmahon(a,b,c), lcm(*range(c+1,c+6)))
    def f(a,b,c):
        return macmahon(a,b,c)/comb(c+a+b-1,a+b-1)
    table = [
        (f(3,b,b),lcm(*range(1,b+3)),f(4,b,b),lcm(*range(1,2*b)))
        for b in (4,5,6)]
    assert table == [(Fraction(588,5),60,Fraction(3528,5),420),
                     (924,420,13068,2520),
                     (Fraction(56628,7),840,Fraction(2024451,7),27720)]
    assert f(3,7,7) == Fraction(306735,4) and f(3,3,6) == 90
    return table


def check_editorial_examples():
    # The fifth-root residual condition and the seventh-root obstruction for the 3-cube.
    assert (q-q**6) - (q**2-q**7) == q*(1-q)*(1-q**5)
    exps = multiplicities(3,3,3)
    assert {d:e for d,e in exps.items() if e} == {2:1,5:1,6:3,7:2,8:1}
    cube_mod_seven = ordinary(3,3,3).applyfunc(lambda x: x % 7)
    assert cube_mod_seven == Matrix([[6,1,6],[1,6,1],[6,1,6]])
    assert integer_factors(Matrix([[2,2],[0,2]])) == (2,2)
    # The two nonconsecutive indices share only Phi_5 and Phi_10 over Q.
    es = elementary(10)
    assert strip_unit(es[1].gcd(es[4])) == phi(5)*phi(10)
    assert gcd(comb(10,1),comb(10,4)) == 10
    assert gcd(*(comb(10,k) for k in range(1,5))) == 5
    # Three active levels can leave two distinct decompositions consistent with the bounds.
    case = (4,13,18)
    exps = multiplicities(*case)
    assert {d:exps[d] for d in (2,4,8,16,32) if exps[d]} == {8:2,16:2,32:3}
    smith = integer_factors(ordinary(*case))
    two_parts = [int(n & -n) for n in smith]
    assert two_parts == [1,2,8,8]
    print("Three-level example",case,"2-primary orders including unit",two_parts)


def main():
    print("Section 6 table", check_examples(), flush=True)
    check_editorial_examples()
    cases=[(1,1,1),(1,4,5),(3,1,2),(3,2,1),(2,2,2),(3,4,3),(4,3,5),(4,4,4),(5,3,3)]
    total=0
    for box in cases:
        n, factors=polynomial_minors(*box)
        print("minors", box, n, "integer Smith factors", factors, flush=True)
        total+=n
    assert total == 453
    print("all polynomial minors", total, flush=True)
    intervals=interval_checks()
    assert intervals == 2275
    print("interval/prime comparisons", intervals, flush=True)
    boxes=list(combinations_with_replacement(range(1,13),3))
    boxes += [(2,1,7),(2,5,22),(2,30,30),(2,8,17),(2,7,15),
              (6,2,4),(9,1,1),(5,3,3),(3,4,3),(4,3,5),(17,17,17)]
    pairs, at_most_two=0,0
    for box in boxes:
        tested, classified=primary_checks(*box)
        pairs+=tested
        at_most_two+=classified
    print("integer Smith forms",len(boxes),"box/prime comparisons",pairs,
          "with at most two active powers",at_most_two)
    print("All checks passed.")


if __name__ == "__main__":
    main()
