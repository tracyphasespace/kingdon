#test_ga_core.py
# tests/test_ga_core.py

import sys, os
import importlib
import pytest

# 1. Import your forked Algebra as GA_FORK
from kingdon.algebra import Algebra as GA_FORK

# 2. Temporarily add the original-kingdon directory to sys.path
orig_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'original_kingdon'))
sys.path.insert(0, orig_dir)
GA_ORIG = importlib.import_module('kingdon.algebra').Algebra
# then put your fork back on top
sys.path.pop(0)

@pytest.fixture(params=[GA_FORK, GA_ORIG], ids=['fork','orig'])
def alg_cls(request):
    return request.param

def test_basis_squares(alg_cls):
    """e_i^2 == signature_i."""
    alg = alg_cls(3,3)
    # signature array: +1 / -1 / 0
    for i, sig in enumerate(alg.signature):
        ei = alg.multivector(name='e', keys=(1<<i,))
        sq = alg.gp(ei, ei)
        # Only scalar part should survive.
        assert set(sq.keys()) <= {0}
        val = sq.values()[0]
        assert pytest.approx(val) == sig

def test_anticommutator(alg_cls):
    """e_i e_j + e_j e_i == 0 for i!=j."""
    alg = alg_cls(3,3)
    for i in range(alg.d):
        for j in range(alg.d):
            if i == j: continue
            ei = alg.multivector(name='e', keys=(1<<i,))
            ej = alg.multivector(name='e', keys=(1<<j,))
            s = alg.add(alg.gp(ei, ej), alg.gp(ej, ei))
            # all components must be zero
            assert all(pytest.approx(0)==v for v in s.values())

def test_outer_product(alg_cls):
    """e_i ^ e_j == e_{ij}."""
    alg = alg_cls(3,3)
    for i in range(alg.d):
        for j in range(i+1, alg.d):
            ei = alg.multivector(name='e', keys=(1<<i,))
            ej = alg.multivector(name='e', keys=(1<<j,))
            op = alg.op(ei, ej)  # exterior product
            # Should produce exactly one key = i|j
            assert op.keys() == ( (1<<i) ^ (1<<j), )
            assert pytest.approx(1) == op.values()[0]

def test_reverse(alg_cls):
    """Reverse(x*y) == reverse(y) * reverse(x)."""
    alg = alg_cls(3,3)
    # pick two random multivectors
    a = alg.multivector(e0=2.0, e12=3.0)
    b = alg.multivector(e1=1.0, e23=4.0)
    lhs = alg.reverse(alg.gp(a, b))
    rhs = alg.gp(alg.reverse(b), alg.reverse(a))
    # keys and values must match
    assert set(lhs.keys()) == set(rhs.keys())
    for k in lhs.keys():
        assert pytest.approx(lhs.items().__dict__['value'] if False else dict(lhs.items())[k]) \
               == dict(rhs.items())[k]

def test_associativity(alg_cls):
    """(a b) c == a (b c)."""
    alg = alg_cls(3,3)
    a = alg.multivector(e0=1, e1=2)
    b = alg.multivector(e2=3, e3=4)
    c = alg.multivector(e12=5, e23=6)
    left = alg.gp(alg.gp(a,b), c)
    right= alg.gp(a, alg.gp(b,c))
    assert set(left.keys()) == set(right.keys())
    for k in left.keys():
        assert pytest.approx(dict(left.items())[k]) == dict(right.items())[k]
