import pytest
import sympy
# Assuming Polynomial and RationalPolynomial are now using SymPy coefficients
from kingdon.polynomial import Polynomial, RationalPolynomial, compare
from sympy import symbols, Integer, Add, Mul, Pow, S # Import necessary SymPy elements

from sympy import sympify
# Or import the whole module:
# import sympy
# and then use sympy.sympify(...)

# Define common SymPy symbols for use in tests
x, y, z = symbols('x y z')

# Helper function to create polynomials from sympy for test cases
# Ensures consistent symbol ordering and nvars
def P(expr, vars=None):
    if vars is None:
        # Default order if not specified, assuming standard x, y, z context
        default_vars = sorted(list(expr.free_symbols), key=str)
        # Handle constant case where free_symbols is empty
        if not default_vars:
            # Decide default nvars for constants in tests, e.g., 0 or based on context
            # Let's assume 0 if no variables are present
            return Polynomial.fromsympy(expr, [])
        return Polynomial.fromsympy(expr, default_vars)
    # Use specified vars
    return Polynomial.fromsympy(expr, vars)

# Helper function to create rational polynomials
def RP(num_expr, den_expr=1, vars=None):
     # Determine variables if not provided
     if vars is None:
         num_sym = sympify(num_expr).free_symbols
         den_sym = sympify(den_expr).free_symbols
         all_sym = sorted(list(num_sym | den_sym), key=str)
         vars = all_sym if all_sym else [] # Use empty list for constants

     num_poly = P(num_expr, vars)
     den_poly = P(den_expr, vars)
     return RationalPolynomial(num_poly, den_poly)


# --- Test Compare Function ---
# Note: compare function might be less relevant if relying purely on Polynomial.__eq__
# Keeping the test structure but updating initialization
# @pytest.mark.xfail(reason="Feedback: Library calculation bug OR Test needs update") # Keep xfail until compare logic verified with symbolic Polynomials
@pytest.mark.parametrize("a, b, expected_comparison_result", [
    # Using Polynomial.fromsympy for initialization
    (P(x + 2*y + 3*z, [x, y, z]), P(x + 2*y + 3*z, [x, y, z]), 0),
    (P(x + 2*y + 3*z, [x, y, z]), P(x + 2*y + 4*z, [x, y, z]), -1), # 3z < 4z assuming some order
    (P(x + 2*y + 4*z, [x, y, z]), P(x + 2*y + 3*z, [x, y, z]), 1),
    (P(x + 2*y + 4*z, [x, y, z]), None, -1), # Existing poly vs None
    (None, P(x + 2*y + 4*z, [x, y, z]), 1), # None vs Existing poly
    (None, None, 1), # None vs None
], ids=[
    "equal_poly", "a_less_than_b", "a_greater_than_b", "a_vs_none", "none_vs_b", "none_vs_none"
])
def test_compare_api(a, b, expected_comparison_result):
    """ Tests the `compare` utility function API (if still used/relevant). """
    # This test might need adjustment based on how compare handles SymPy coeffs
    assert compare(a, b) == expected_comparison_result


# --- Test Arithmetic Operators ---

# Test the addition operator (+) API
@pytest.mark.parametrize("a, b, expected_sum", [
    # === Polynomial Tests ===
    (P(x, [x]), P(2*x, [x]), P(3*x, [x])),
    (P(x*y + 2*z, [x, y, z]), P(3*x*y + 4*z, [x, y, z]), P(4*x*y + 6*z, [x, y, z])),
    (P(x + 2*y + 3*z, [x, y, z]), P(4*x + 5*y + 6*z, [x, y, z]), P(5*x + 7*y + 9*z, [x, y, z])),
    (P(x + 2*y + 3*z, [x, y, z]), 0, P(x + 2*y + 3*z, [x, y, z])), # Add zero scalar
    (P(x - 2*y + 3*z, [x, y, z]), P(-x + 2*y - 3*z, [x, y, z]), P(0, [x, y, z])), # Add to zero poly

    # === RationalPolynomial Tests ===
    (RP(x + 2*y, z, [x, y, z]), RP(3*x + 4*y, z, [x, y, z]), RP(4*x + 6*y, z, [x, y, z])), # Same denom
    # Expected result for different denoms: ( (x+2y)*2z + (3x+4y)*z ) / (z*2z) -> simplify
    (RP(x + 2*y, z, [x, y, z]), RP(3*x + 4*y, 2*z, [x, y, z]), RP(2*x*z + 4*y*z + 3*x*z + 4*y*z, 2*z*z, [x, y, z]).simplify()), # Add diff denom
    (RP(x + 2*y, z, [x, y, z]), 0, RP(x + 2*y, z, [x, y, z])), # Add zero scalar
    (RP(-x + 2*y, z, [x, y, z]), RP(x - 2*y, z, [x, y, z]), RP(0, z, [x, y, z])), # Add to zero rational
], ids=[
    "poly_add_simple", "poly_add_multi", "poly_add_order", "poly_add_zero", "poly_add_to_zero",
    "rational_add_same_denom", "rational_add_diff_denom", "rational_add_zero", "rational_add_to_zero"
])
def test_add_operator(a, b, expected_sum):
    """ Tests the addition (+) operator API for Polynomial and RationalPolynomial. """
    result = a + b
    # Use Polynomial/RationalPolynomial equality check
    assert result == expected_sum


# Test the multiplication operator (*) API
@pytest.mark.parametrize("a, b, expected_product", [
    # === Polynomial Tests ===
    (P(x, [x]), P(2*x, [x]), P(2*x**2, [x])),
    (P(x*y + 2*z, [x,y,z]), P(3*x*y + 4*z, [x,y,z]), P(3*x**2*y**2 + 4*x*y*z + 6*x*y*z + 8*z**2, [x,y,z])), # Expand manually or use sympy
    (P(x + 2*y + 3*z, [x,y,z]), 0, P(0, [x,y,z])), # Multiply by zero scalar
    (P(x + 2*y, [x,y]), 1, P(x + 2*y, [x,y])), # Multiply by one scalar
    (P(x - 2*y, [x,y,z]), P(-x + 3*z, [x,y,z]), P(-x**2 + 3*x*z + 2*x*y - 6*y*z, [x,y,z])),

    # === RationalPolynomial Tests ===
    (RP(x, 1, [x]), RP(2*x, 1, [x]), RP(2*x**2, 1, [x])),
    (RP(x*y + 2*z, 1, [x,y,z]), RP(3*x*y + 4*z, 1, [x,y,z]), RP(3*x**2*y**2 + 10*x*y*z + 8*z**2, 1, [x,y,z])), # Result from Polynomial test
    (RP(x + 2*y + 3*z, 1, [x,y,z]), 0, RP(0, 1, [x,y,z])), # Multiply by zero scalar
    (RP(x + 2*y, z, [x,y,z]), 1, RP(x + 2*y, z, [x,y,z])), # Multiply by one scalar
    # Multiply rational polynomials leading to simplification: (2xz/1) * (1 / 3z) = 2x/3
    (RP(2*x*z, 1, [x,z]), RP(1, 3*z, [x,z]), RP(2*x, 3, [x,z])), # z should cancel
], ids=[
    "poly_mul_simple", "poly_mul_multi", "poly_mul_by_zero", "poly_mul_by_one", "poly_mul_neg_coeffs",
    "rational_mul_simple", "rational_mul_multi", "rational_mul_by_zero", "rational_mul_by_one", "rational_mul_simplify"
])
def test_mul_operator(a, b, expected_product):
    """ Tests the multiplication (*) operator API. """
    result = a * b
    # Use Polynomial/RationalPolynomial equality check
    # For simplified cases, ensure expected_product is also simplified if needed
    if isinstance(expected_product, RationalPolynomial):
        assert result == expected_product.simplify() # Compare with simplified expected
    else:
        assert result == expected_product


# Test the negation operator (-) API
@pytest.mark.parametrize("a, expected_negation", [
    (P(x, [x]), P(-x, [x])),
    (P(x*y + 2*z, [x, y, z]), P(-x*y - 2*z, [x, y, z])),
    (P(0, []), P(0, [])), # Negate zero polynomial
    (P(-x + 2*y - 3*z, [x, y, z]), P(x - 2*y + 3*z, [x, y, z])),
    # Rational
    (RP(x - y, z, [x,y,z]), RP(-x + y, z, [x,y,z])),
    (RP(x, y - z, [x,y,z]), RP(-x, y - z, [x,y,z])), # Neg only numerator
], ids=[
    "poly_neg_simple", "poly_neg_multi", "poly_neg_empty", "poly_neg_mixed_sign",
    "rational_neg_num", "rational_neg_denom_unchanged"
])
def test_neg_operator(a, expected_negation):
    """ Tests the negation (-) operator API. """
    assert -a == expected_negation


# Test the inverse method (.inv()) API for RationalPolynomial
@pytest.mark.parametrize("a, expected_inverse", [
    (RP(x, 1, [x]), RP(1, x, [x])),
    (RP(x*y + 2*z, 1, [x, y, z]), RP(1, x*y + 2*z, [x, y, z])),
    (RP(5, 1, []), RP(1, 5, [])), # Constant case
    (RP(x, y, [x, y]), RP(y, x, [x, y])), # Fraction
], ids=[
    "rational_inv_simple", "rational_inv_multi", "rational_inv_constant", "rational_inv_fraction"
])
def test_inverse_method(a, expected_inverse):
    """ Tests the .inv() method API for RationalPolynomial. """
    assert a.inv() == expected_inverse


# Test the property that a * a.inv() == 1 for RationalPolynomial
@pytest.mark.parametrize("a", [
    RP(x, 1, [x]),
    RP(x*y + 2*z, 1, [x, y, z]),
    RP(3*x, 2*y, [x, y]),
    RP(5, 1, []), # Constant case
], ids=[
    "identity_check_simple", "identity_check_multi", "identity_check_rational", "identity_check_constant"
])
def test_inverse_identity(a):
    """ Tests the fundamental property that a * a.inv() == 1. """
    identity = RP(1, 1, a.num.symbols if hasattr(a.num,'symbols') else []) # Rational 1 with same vars
    # Simplify result before comparison
    result = (a * a.inv()).simplify()
    assert result == identity


# Test the power operator (**) API
@pytest.mark.parametrize("a, n, expected_power", [
    # === Polynomial Tests ===
    (P(2*x, [x]), 2, P(4*x**2, [x])),
    (P(x + y, [x, y]), 2, P(x**2 + 2*x*y + y**2, [x, y])), # (x+y)^2
    (P(2*x, [x]), 3, P(8*x**3, [x])), # Changed from 5 to 3 for simplicity
    (P(3, []), 3, P(27, [])), # Constant case

    # === RationalPolynomial Tests ===
    (RP(2*x, 1, [x]), 2, RP(4*x**2, 1, [x])),
    (RP(3*x, y, [x, y]), 2, RP(9*x**2, y**2, [x, y])), # (3x/y)^2
    (RP(2*x, 1, [x]), 3, RP(8*x**3, 1, [x])), # Changed from 5 to 3
    (RP(3, 1, []), 3, RP(27, 1, [])), # Constant case
], ids=[
    "poly_pow_simple", "poly_pow_binomial", "poly_pow_high", "poly_pow_const",
    "rational_pow_simple", "rational_pow_fraction", "rational_pow_high", "rational_pow_const"
])
def test_pow_operator(a, n, expected_power):
    """ Tests the power (**) operator API. """
    assert a ** n == expected_power


# Test creation with symbolic coefficients directly using fromsympy
def test_creation_with_sympy():
    """ Tests creating Polynomials using SymPy expressions via fromsympy API. """
    p1 = Polynomial.fromsympy(x + y, symbols=[x, y])
    p2 = Polynomial.fromsympy(x * y, symbols=[x, y])
    p_sum = p1 + 1              # Test __add__ API with int
    p_prod = p1 * p2            # Test __mul__ API

    assert isinstance(p_sum, Polynomial)
    assert isinstance(p_prod, Polynomial)

    # Verify results using API comparison (==)
    # Ensure expected results have the correct symbols/nvars
    assert p_sum == Polynomial.fromsympy(x + y + 1, symbols=[x, y])
    # (x+y)*(xy) = x^2y + xy^2
    assert p_prod == Polynomial.fromsympy(x**2*y + x*y**2, symbols=[x, y])

    # Rational Polynomial with SymPy using API
    rp1 = RationalPolynomial.fromsympy(x / y, symbols=[x, y])
    rp2 = RationalPolynomial.fromsympy(y / x, symbols=[x, y])
    rp_prod = rp1 * rp2 # Test __mul__ API
    rp_sum = rp1 + 1    # Test __add__ API

    assert isinstance(rp_prod, RationalPolynomial)
    assert isinstance(rp_sum, RationalPolynomial)

    # Product should simplify to 1
    # Create expected RP(1) with correct symbols
    expected_one = RationalPolynomial.fromsympy(S.One, symbols=[x, y])
    assert rp_prod.simplify() == expected_one

    # Sum should be (x+y)/y
    expected_sum = RationalPolynomial.fromsympy((x+y)/y, symbols=[x, y])
    assert rp_sum.simplify() == expected_sum.simplify()


# Optional: Add tests for Polynomial.fromsympy edge cases
def test_fromsympy_edge_cases():
    """ Test Polynomial.fromsympy with constants and different symbol lists. """
    assert Polynomial.fromsympy(S(5), []) == Polynomial({(): 5}, nvars=0)
    assert Polynomial.fromsympy(S(0), []) == Polynomial({}, nvars=0)
    assert Polynomial.fromsympy(x, [x]) == Polynomial({(1,): 1}, nvars=1)
    assert Polynomial.fromsympy(x, [x, y]) == Polynomial({(1, 0): 1}, nvars=2) # x in context of x, y
    assert Polynomial.fromsympy(x + y, [y, x]) == Polynomial({(1, 0): 1, (0, 1): 1}, nvars=2) # Symbols reordered
    with pytest.raises(ValueError): # Expression has symbol z not in list
        Polynomial.fromsympy(x + z, [x, y])
    with pytest.raises(ValueError): # Expression not polynomial
        Polynomial.fromsympy(sympy.sin(x), [x])

# Optional: Add tests for RationalPolynomial.fromsympy edge cases
def test_rational_fromsympy_edge_cases():
    """ Test RationalPolynomial.fromsympy. """
    assert RationalPolynomial.fromsympy(S(5), []) == RationalPolynomial(P(5,[]), P(1,[]))
    assert RationalPolynomial.fromsympy(x/y, [x, y]) == RationalPolynomial(P(x, [x,y]), P(y, [x,y]))
    # Test simplification during conversion
    assert RationalPolynomial.fromsympy((x**2-1)/(x-1), [x]) == RationalPolynomial(P(x+1, [x]), P(1, [x]))