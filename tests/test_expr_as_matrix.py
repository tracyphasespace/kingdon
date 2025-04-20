import pytest
from sympy import Matrix, symbols
import numpy as np
import sympy

from kingdon import Algebra, MultiVector
# Assuming expr_as_matrix is part of the intended public API for matrix representations
from kingdon.matrixreps import expr_as_matrix

@pytest.fixture(scope='module')
def alg3d():
    """Provides a 3D Euclidean Algebra R(3,0,0) instance."""
    return Algebra(3, 0)

@pytest.fixture
def symbolic_vec(alg3d):
    """Provides a symbolic vector in R(3)."""
    return alg3d.vector(name='x') # Creates vector with symbolic components x_e1, x_e2, x_e3

@pytest.fixture
def symbolic_bivec(alg3d):
    """Provides a symbolic bivector in R(3)."""
    return alg3d.bivector(name='B') # Creates bivector B_e12*e12 + B_e13*e13 + B_e23*e23

@pytest.mark.xfail(reason="Feedback: Feature not implemented")
def test_expr_as_matrix_symbolic_cp(alg3d, symbolic_vec, symbolic_bivec):
    """
    Tests expr_as_matrix API for the commutator product (cp) with symbolic inputs.
    Verifies the matrix for B cp x, where x is a vector.
    """
    B = symbolic_bivec
    x = symbolic_vec
    B12, B13, B23 = B.values() # Extract symbolic coefficients via API

    # --- Test Case 1: Full vector output ---
    # Get matrix A such that A @ x_coeffs = (B cp x)_coeffs, using expr_as_matrix API
    A, y = expr_as_matrix(alg3d.cp, B, x)

    # Verify the resulting matrix A against the expected symbolic matrix
    expected_A = Matrix([[0, B12, B13], [-B12, 0, B23], [-B13, -B23, 0]])
    assert A == expected_A, "Matrix for B cp x (vector output) is incorrect."

    # Verify the resulting multivector y is correct using API comparison
    assert y == B.cp(x), "Resulting multivector y does not match direct calculation."
    assert y.grades == (1,), f"Resulting multivector y should be grade 1, got {y.grades}" # Check grade via API

    # --- Test Case 2: Filtered vector output using res_like API ---
    # Define a target shape for the output (only e1 and e3 components) via API
    X_target = alg3d.vector(e1=1, e3=1) # Example MV with desired output structure
    # Call expr_as_matrix API with res_like argument
    A_filtered, y_filtered = expr_as_matrix(alg3d.cp, B, x, res_like=X_target)

    # Verify the filtered matrix A_filtered
    # Should only contain rows corresponding to e1 and e3 outputs
    expected_A_filtered = Matrix([[0, B12, B13], [-B13, -B23, 0]])
    assert A_filtered == expected_A_filtered, "Filtered matrix for B cp x (e1, e3 output) is incorrect."

    # Verify the filtered multivector y_filtered using API checks
    expected_y_filtered_keys = X_target.keys() # Get keys via API
    assert y_filtered.keys() == expected_y_filtered_keys, "Filtered result y_filtered has incorrect keys."
    # Check components match the full result's components for those keys using API attribute access
    assert y_filtered.e1 == y.e1
    assert y_filtered.e3 == y.e3
    assert y_filtered.e2 == 0 # Ensure other components are zero


@pytest.mark.xfail(reason="Feedback: Feature not implemented")
def test_expr_as_matrix_symbolic_acp(alg3d, symbolic_vec, symbolic_bivec):
    """
    Tests expr_as_matrix API for the anti-commutator product (acp) with symbolic inputs.
    Verifies the matrix for B acp x, where x is vector. Expects trivector output in R(3).
    """
    B = symbolic_bivec
    x = symbolic_vec
    B12, B13, B23 = B.values() # Extract symbolic coefficients via API

    # Get matrix A such that A @ x_coeffs = (B acp x)_coeffs, using expr_as_matrix API
    # Input x is grade 1, B is grade 2. B acp x is expected grade 3.
    A, y = expr_as_matrix(alg3d.acp, B, x)

    # Verify the resulting matrix A against the expected symbolic matrix
    # Output is grade 3 (e123), input is grade 1 (e1, e2, e3)
    # Matrix maps 3 inputs to 1 output. Should be 1x3.
    expected_A = Matrix([[B23, -B13, B12]]) # Based on B acp x formula
    assert A == expected_A, "Matrix for B acp x (trivector output) is incorrect."

    # Verify the resulting multivector y is correct using API checks
    assert y == B.acp(x), "Resulting multivector y does not match direct calculation."
    assert y.grades == (3,), f"Resulting multivector y should be grade 3, got {y.grades}"


@pytest.mark.xfail(reason="Feedback: Feature not implemented")
def test_expr_as_matrix_symbolic_sw(alg3d, symbolic_bivec):
    """
    Tests expr_as_matrix API for the sandwich product (sw) with symbolic inputs.
    Verifies the matrix for B sw x (conjugation of vector x BY bivector B).
    """
    B = symbolic_bivec
    x_sym = alg3d.vector(name='x') # Create symbolic vector via API

    # Define the operation (B sw V) as a lambda to pass to expr_as_matrix API
    sandwich_op = lambda v: B.sw(v)
    # Get matrix A such that A @ x_coeffs = (B sw x)_coeffs using expr_as_matrix API
    A, y = expr_as_matrix(sandwich_op, x_sym) # Pass op and the vector being transformed

    # --- Verification ---
    # Verify the resulting multivector y against direct calculation using API
    assert y == B.sw(x_sym), "Resulting multivector y does not match direct calculation."
    assert y.grades == (1,), f"Resulting multivector y should be grade 1 (vector), got {y.grades}" # Check grade via API

    # Verify the matrix A by applying it to a basis vector (e.g., e3) and comparing
    e3 = alg3d.blades.e3
    y_e3_direct = B.sw(e3) # Calculate B sw e3 directly via API

    # Calculate A @ e3_coeffs (coeffs are [0, 0, 1])
    e3_coeffs = Matrix([0, 0, 1])
    A_times_e3_coeffs = A * e3_coeffs

    # Construct the expected result vector from A_times_e3_coeffs
    y_e3_from_matrix = alg3d.vector([A_times_e3_coeffs[0], A_times_e3_coeffs[1], A_times_e3_coeffs[2]])

    # Compare the two methods of getting the result for e3 using API equality
    # Use simplify for robust comparison of symbolic expressions
    diff = y_e3_direct - y_e3_from_matrix
    assert all(sympy.simplify(c) == 0 for c in diff.values()), \
           f"Matrix applied to e3 coeffs does not match direct B sw e3. Matrix result: {y_e3_from_matrix}, Direct: {y_e3_direct}"


@pytest.mark.xfail(reason="Feedback: Feature not implemented")
def test_expr_as_matrix_numerical(alg3d):
    """
    Tests expr_as_matrix API with numerical inputs.
    """
    # Use a specific numerical bivector B = e12 + 2*e13 using API
    B = alg3d.blades.e12 + 2*alg3d.blades.e13
    # Use a specific numerical vector x = e1 + 3*e2 using API
    x = alg3d.blades.e1 + 3*alg3d.blades.e2

    # --- Test Commutator Product using API ---
    A_cp, y_cp = expr_as_matrix(alg3d.cp, B, x)

    # Check matrix type and shape
    assert isinstance(A_cp, np.ndarray), "Numerical matrix A should be np.ndarray"
    # Input x is vector (3 components), output y is vector (3 components). Matrix should be 3x3.
    assert A_cp.shape == (3, 3), f"Matrix A shape incorrect, expected (3, 3), got {A_cp.shape}"

    # Verify the matrix content for this specific B = 1*e12 + 2*e13 + 0*e23
    # B12=1, B13=2, B23=0
    # expected_A = [[0, B12, B13], [-B12, 0, B23], [-B13, -B23, 0]]
    expected_A_cp = np.array([
        [0, 1, 2],
        [-1, 0, 0],
        [-2, 0, 0]
    ])
    assert np.allclose(A_cp, expected_A_cp), "Numerical matrix A content incorrect for B cp x."

    # Verify the resulting multivector y using API equality
    y_direct_cp = B.cp(x) # Direct calculation via API
    assert y_cp == y_direct_cp, "Numerical result y incorrect for B cp x."

    # Verify matrix multiplication relation: A @ x_coeffs = y_coeffs using API components
    x_coeffs = np.array([x.e1, x.e2, x.e3]) # [1, 3, 0]
    y_coeffs = np.array([y_cp.e1, y_cp.e2, y_cp.e3])
    assert np.allclose(A_cp @ x_coeffs, y_coeffs), "Matrix multiplication A @ x != y"


@pytest.mark.xfail(reason="Feedback: Feature not implemented")
def test_expr_as_matrix_numerical_array_input(alg3d):
    """
    Tests expr_as_matrix API with numerical, array-valued inputs.
    """
    # Create an array-valued Bivector (e.g., 2 different bivectors) using API
    B_vals = np.array([
        [1, 2, 0],  # Bivector 1: 1*e12 + 2*e13 + 0*e23
        [0, 1, 3]   # Bivector 2: 0*e12 + 1*e13 + 3*e23
    ]).T # Shape (3, 2) -> coeffs for e12, e13, e23 for 2 MVs
    B_array = alg3d.bivector(B_vals) # Create via API
    assert B_array.shape == (2,) # Check shape via API

    # Create a symbolic vector x (matrix result will depend on x) using API
    x = alg3d.vector(name='x')

    # Get the matrix representation for the commutator product B cp x using API
    # Expect a result where A is a list or array of matrices, one for each MV in B_array
    A_list, y_array = expr_as_matrix(alg3d.cp, B_array, x)

    # Check the type of A (should be list or array of symbolic matrices)
    assert isinstance(A_list, list) or isinstance(A_list, np.ndarray), \
           "Matrix result A should be list or ndarray for array input."
    assert len(A_list) == 2, "Should get 2 matrices for 2 input bivectors."

    # Check the type of the resulting multivector (should be array-valued)
    assert isinstance(y_array, MultiVector)
    assert y_array.shape == (2,) # Check shape via API

    # Verify the first matrix corresponds to the first bivector (using API to get values)
    B1_coeffs = B_array[0].values() # Values for the first bivector via API __getitem__
    B12_1, B13_1, B23_1 = B1_coeffs[0], B1_coeffs[1], B1_coeffs[2] # Coeffs for e12, e13, e23
    expected_A1 = Matrix([[0, B12_1, B13_1], [-B12_1, 0, B23_1], [-B13_1, -B23_1, 0]])
    assert A_list[0] == expected_A1

    # Verify the second matrix corresponds to the second bivector (using API to get values)
    B2_coeffs = B_array[1].values() # Values for the second bivector via API __getitem__
    B12_2, B13_2, B23_2 = B2_coeffs[0], B2_coeffs[1], B2_coeffs[2] # Coeffs for e12, e13, e23
    expected_A2 = Matrix([[0, B12_2, B13_2], [-B12_2, 0, B23_2], [-B13_2, -B23_2, 0]])
    assert A_list[1] == expected_A2

    # Verify the resulting array multivector using API equality
    assert y_array == B_array.cp(x), "Resulting array multivector y does not match direct calculation."