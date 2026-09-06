#!/usr/bin/env python3
"""Rank-one supported sheaf indices and magnetic-probe pairings for X9.

Exact standard-library checks. This computes the Hilbert-scheme component
of the compact pure-sheaf moduli at fixed surface line bundle, not the full
physical BPS index at the cooperative singularity. The geometric proof of
the component identification is in the manuscript, not established by this
finite arithmetic check.

Run: python3 checks/check_x9_rank_one_index.py
"""

import contextlib
import io
from collections import defaultdict
from fractions import Fraction as Fr
from itertools import product

with contextlib.redirect_stdout(io.StringIO()):
    import check_x9_smooth_side as bulk
    import check_x9_m5_string_anomaly as anomaly

ORDER = 12


def multiply(a, b, order=ORDER):
    result = defaultdict(int)
    for (n, j), c in a.items():
        for (m, k), d in b.items():
            if n + m <= order:
                result[n + m, j + k] += c * d
    return dict(result)


def oscillator_product(b2):
    """Centered Poincare product, exponents recording q and y."""
    series = {(0, 0): 1}
    for m in range(1, ORDER + 1):
        for weight in [-2] + [0] * b2 + [2]:
            factor = {(m * k, weight * k): 1 for k in range(ORDER // m + 1)}
            series = multiply(series, factor)
    return series


def partitions_by_enumeration(n, largest=None):
    """Enumerate ordinary partitions independently of the oscillator product."""
    if n == 0:
        yield ()
        return
    if largest is None:
        largest = n
    for k in range(min(n, largest), 0, -1):
        for rest in partitions_by_enumeration(n - k, k):
            yield (k,) + rest


partition_numbers = [sum(1 for _ in partitions_by_enumeration(n))
                     for n in range(ORDER + 1)]


def colored_partitions(colors):
    series = [1] + [0] * ORDER
    for _ in range(colors):
        series = [sum(series[k] * partition_numbers[n-k] for k in range(n+1))
                  for n in range(ORDER+1)]
    return series


before, after = oscillator_product(2), oscillator_product(1)
exceptional = {(n, 0): p for n, p in enumerate(partition_numbers)}
assert before == multiply(after, exceptional)

counts = {}
for name, series, chi in [('F1', before, 4), ('P2', after, 3)]:
    counts[name] = [sum(c for (m, j), c in series.items() if m == n)
                    for n in range(ORDER + 1)]
    assert counts[name] == colored_partitions(chi)
    for (n, j), coefficient in series.items():
        assert coefficient > 0 and j % 2 == 0 and abs(j) <= 2*n
        assert series.get((n, -j), 0) == coefficient
    b2 = chi - 2
    assert {j: c for (n, j), c in series.items() if n == 1} == {-2: 1, 0: b2, 2: 1}
    assert {j: c for (n, j), c in series.items() if n == 2} == {
        -4: 1, -2: b2+1, 0: 1+b2+b2*(b2+1)//2, 2: b2+1, 4: 1}

assert counts['F1'][:7] == [1, 4, 14, 40, 105, 252, 574]
assert counts['P2'][:7] == [1, 3, 9, 22, 51, 108, 221]
assert Fr(-4, 24) - Fr(-3, 24) == -Fr(1, 24)  # eta^{-1} vacuum shift


def dot(a, b, signature):
    return sum(x*y*s for x, y, s in zip(a, b, signature))


# GRR: Gamma=ch(i_*F)sqrt(Td(X)); F=I_Z tensor L, l=c1(L).
# Q2=i_*(l-K/2), Q0=(l-K/2)^2/2+chi(S)/24-n.
for signature, canonical, chi, restrictions in [
        ((1, -1), (-3, 1), 4, anomaly.restr_before),
        ((1,), (-3,), 3, anomaly.restr_after)]:
    k2 = dot(canonical, canonical, signature)
    c2x = chi-k2
    for coeff in product(range(-3, 4), repeat=len(signature)):
        flux = tuple(Fr(x)-Fr(k, 2) for x, k in zip(coeff, canonical))
        for n in range(4):
            q0 = dot(flux, flux, signature)/2 + Fr(chi, 24)-n
            grr = (Fr(dot(coeff, coeff, signature), 2)
                   - Fr(dot(coeff, canonical, signature), 2)
                   + Fr(k2, 6) + Fr(c2x, 24)-n)
            assert q0 == grr
        electric = [dot(d, flux, signature) for d in restrictions]
        u = Fr(coeff[0])+Fr(3, 2)
        expected = [0, u, -3*u]
        if len(signature) == 2:
            v = Fr(coeff[1])-Fr(1, 2)
            expected[2] -= v
            assert v.denominator == 2 and v != 0
            assert electric[2].denominator == 1
        else:
            assert electric[2].denominator == 2
        assert electric == expected

# A rank-one D4 core has no D6 charge. Against (0,0,gamma,n), its
# Dirac pairing is S.gamma, independent of induced D2,D0 flux charges.
curve_columns = [
    [bulk.tri(d, bulk.D6, bulk.D7) for d in
     [bulk.D2, bulk.D5, bulk.D6, bulk.D7, bulk.D8, bulk.E]],
    [bulk.tri(d, bulk.D2, bulk.D4) for d in
     [bulk.D2, bulk.D5, bulk.D6, bulk.D7, bulk.D8, bulk.E]]]
root = [a+b for a, b in zip(*curve_columns)]
probes = {'D8': bulk.D8, 'D5': bulk.D5, 'D6': bulk.D6}
probe_pairings = {name: [sum(x*y for x, y in zip(divisor, c))
                         for c in curve_columns+[root]]
                  for name, divisor in probes.items()}
assert probe_pairings == {'D8': [0, 0, 0], 'D5': [0, 1, 1], 'D6': [-1, 0, -1]}
for a, b in product(range(-6, 7), repeat=2):
    charge = [a*x+b*y for x, y in zip(*curve_columns)]
    assert sum(x*y for x, y in zip(bulk.D8, charge)) == 0
for d in (bulk.D5, bulk.D6):
    # These restrictions vanish on the full F1 Picard lattice, generated by
    # h=D1|S and e=E|S; hence every supported flux has both effective charges zero.
    assert bulk.tri(bulk.D8, d, bulk.D1) == 0
    assert bulk.tri(bulk.D8, d, bulk.E) == 0

print('Fixed-line-bundle component: Hilb^n(S), numerical DT sign (+1).')
print('F1 Euler indices n=0,...,12:', counts['F1'])
print('P2 Euler indices n=0,...,12:', counts['P2'])
print('Refined ratio through q^12: product_m (1-q^m)^(-1), independent of y.')
print('Vacuum-normalized unrefined factors: eta^(-4), eta^(-3).')
print('GRR and half-canonical flux shifts checked in both surface lattices.')
print('Probe pairings with (C1,C2,R):', probe_pairings)
print('D8 ordinary halo wall-crossing multiplier from these charges is trivial.')
print('ALL CHECKS PASSED; no full compact-condensate index or threshold no-go asserted.')
