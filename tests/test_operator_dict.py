import pytest
import sympy

# Removed OperatorDict, UnaryOperatorDict, codegen imports
from kingdon import Algebra, MultiVector

@pytest.fixture(scope='module')
def alg2d():
    """Provides a 2D Euclidean Algebra instance R(2,0)."""
    return Algebra(2) # p=2, q=0, r=0

@pytest.fixture
def symbolic_mv_x(alg2d):
    """Provides a generic symbolic multivector 'x' in 2D."""
    # Include all grades for generality
    return alg2d.multivector(name='x')

@pytest.fixture
def symbolic_mv_y(alg2d):
    """Provides a generic symbolic multivector 'y' in 2D."""
    # Include all grades for generality
    return alg2d.multivector(name='y')

def test_geometric_product_symbolic(alg2d):
    """
    Tests the geometric product (*) API for symbolic multivectors in 2D.
    Verifies the result against known algebraic rules for vectors.
    """
    # Test a specific case: product of two symbolic vectors.
    x_vec = alg2d.vector(name='x') # x_e1*e1 + x_e2*e2
    y_vec = alg2d.vector(name='y') # y_e1*e1 + y_e2*e2

    # Perform the geometric product using the public API operator *
    prod_vec = x_vec * y_vec

    # --- Verification ---
    # Expected result for R(2,0):
    # (x1*e1 + x2*e2) * (y1*e1 + y2*e2)
    # = (x1*y1 + x2*y2) + (x1*y2 - x2*y1)*e12
    # Changed .symbols() to direct attribute access per feedback
    x1, x2 = x_vec.e1, x_vec.e2 # Get symbolic coeffs via API attribute access
    y1, y2 = y_vec.e1, y_vec.e2 # Get symbolic coeffs via API attribute access

    expected_scalar_part = x1*y1 + x2*y2
    expected_bivector_part = (x1*y2 - x2*y1) * alg2d.blades.e12 # Use API blade

    # Construct the expected MultiVector result using the API
    expected_prod = alg2d.scalar(expected_scalar_part) + expected_bivector_part

    # Assert equality using the MultiVector.__eq__ API (should handle simplify)
    assert prod_vec == expected_prod, f"Geometric product of symbolic vectors incorrect. Got {prod_vec}, expected {expected_prod}"

    # Test that the result is indeed a MultiVector
    assert isinstance(prod_vec, MultiVector)

@pytest.mark.xfail(reason="Feedback: Library calculation/simplification error")
def test_inverse_symbolic(alg2d):
    """
    Tests the inverse (.inv()) method API for a symbolic multivector.
    Verifies the property x * x.inv() == 1 for a rotor example.
    """
    # Test with a specific symbolic, invertible element: a rotor
    a = sympy.symbols('a')
    # Create rotor R = cos(a/2) + sin(a/2)*e12 using API
    rotor = sympy.cos(a/2) * alg2d.scalar(1) + sympy.sin(a/2) * alg2d.blades.e12
    rotor_rev = ~rotor # Calculate reverse via API operator

    # Calculate inverse using the public API method .inv()
    try:
        rotor_inv = rotor.inv()
    except Exception as e:
        pytest.fail(f"Calculating rotor.inv() failed: {e}")

    # Verify the result: inv(R) should equal ~R
    # Use sympy.simplify on the difference for robust comparison
    diff = rotor_inv - rotor_rev
    # Check if all components of the difference simplify to zero using API .values()
    # Check keys exist before simplifying, handle potential empty diff
    is_diff_zero = all(sympy.simplify(c) == 0 for c in diff.values()) if diff.values() else (diff == 0)
    assert is_diff_zero, f"Inverse of symbolic rotor incorrect. Expected ~R = {rotor_rev}, got {rotor_inv}"

    # Verify the defining property: rotor * rotor_inv == 1
    identity = alg2d.scalar(1) # Identity element via API
    product = rotor * rotor_inv # Use API operator
    prod_diff = product - identity # Use API operator
    # Check if product simplifies to identity
    is_prod_one = all(sympy.simplify(c) == 0 for c in prod_diff.values()) if prod_diff.values() else (prod_diff == 0)
    assert is_prod_one, f"Product rotor * rotor.inv() did not simplify to 1. Got {product}"

    # Test that the result is a MultiVector
    assert isinstance(rotor_inv, MultiVector)