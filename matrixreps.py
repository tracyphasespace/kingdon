"""
Matrix Representations Module for Geometric Algebra
==================================================

This module provides functions to create matrix representations of Geometric Algebra
basis blades. These matrix representations allow Geometric Algebra operations to be
performed using standard matrix operations.

The primary function is `matrix_rep()`, which constructs matrix representations
for a given signature (p,q,r).

Functions:
    matrix_rep: Creates matrix representations of basis blades for a Clifford algebra
    kronecker_product: Computes the Kronecker product of matrices
    matrix_rep_blade: Creates the matrix representation for a specific blade

Example:
    >>> from kingdon.matrixreps import matrix_rep
    >>> # Create matrix representations for a 3D Euclidean geometric algebra
    >>> matrices = matrix_rep(p=3)
    >>> # matrices[0] is the identity (scalar), matrices[1] is e1, matrices[2] is e2, etc.
"""

from typing import List, Tuple, Optional, Union, Sequence, Dict, Any, cast
import numpy as np
from functools import reduce
import operator

def expr_as_matrix(op, mv, res_like=None):
    raise NotImplementedError("expr_as_matrix not implemented in this fork")


def kronecker_product(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """
    Compute the Kronecker product of two matrices.
    
    The Kronecker product is a generalization of the outer product from vectors to matrices.
    
    Args:
        A: First input matrix
        B: Second input matrix
        
    Returns:
        The Kronecker product A ⊗ B
        
    Example:
        >>> A = np.array([[1, 2], [3, 4]])
        >>> B = np.array([[0, 1], [1, 0]])
        >>> kronecker_product(A, B)
        array([[0, 1, 0, 2],
               [1, 0, 2, 0],
               [0, 3, 0, 4],
               [3, 0, 4, 0]])
    """
    return np.kron(A, B)

def matrix_rep(p: int = 0, q: int = 0, r: int = 0,
               custom_matrices: Optional[Dict[str, np.ndarray]] = None) -> List[np.ndarray]:
    """
    Create matrix representations for the basis blades of a Clifford algebra.

    This implementation uses Kronecker products of 2x2 matrices to build representations
    for each basis blade. The resulting matrices can be used to perform Geometric Algebra
    operations using standard matrix operations.
    
    Args:
        p: Number of positive dimensions (default: 0)
        q: Number of negative dimensions (default: 0)
        r: Number of null dimensions (default: 0)
        custom_matrices: Optional dictionary of custom matrices to use for 
                         positive, negative, and null dimensions
    
    Returns:
        A list of numpy arrays representing each basis blade.
        The list is ordered by binary count, with the scalar (1) first.
    
    Example:
        >>> mats = matrix_rep(p=3, q=0, r=0)
        >>> mats[0]  # Identity matrix (scalar)
        >>> mats[1]  # First basis vector e₁
        >>> mats[3]  # Bivector e₁₂
        
    Notes:
        The total number of matrices returned is 2^d where d = p + q + r.
        The matrices use the standard 2x2 matrix representations:
        - Positive dimensions: [[0, 1], [1, 0]]
        - Negative dimensions: [[0, 1], [-1, 0]]
        - Null dimensions: [[0, 0], [1, 0]]
    """
    d = p + q + r
    
    # Define standard 2x2 matrices
    I2 = np.array([[1, 0], [0, 1]])  # Identity
    P2 = np.array([[0, 1], [1, 0]])  # Positive signature (squares to +1)
    N2 = np.array([[0, 1], [-1, 0]]) # Negative signature (squares to -1)
    Z2 = np.array([[0, 0], [1, 0]])  # Null dimensions (squares to 0)
    
    # Allow for custom matrices if provided
    if custom_matrices:
        P2 = custom_matrices.get('positive', P2)
        N2 = custom_matrices.get('negative', N2)
        Z2 = custom_matrices.get('null', Z2)

    # Build list of matrices according to signature
    signature_matrices = [Z2] * r + [P2] * p + [N2] * q
    
    # Create representations for all basis blades
    basis_matrices: List[np.ndarray] = []
    for i in range(2 ** d):
        basis_matrices.append(matrix_rep_blade(i, d, signature_matrices, I2))
    
    return basis_matrices

def matrix_rep_blade(blade_index: int, dimension: int, 
                     signature_matrices: List[np.ndarray], 
                     identity: np.ndarray) -> np.ndarray:
    """
    Create the matrix representation for a specific basis blade.
    
    Args:
        blade_index: Binary index representing the blade
        dimension: Total dimension of the algebra (p+q+r)
        signature_matrices: List of matrices representing each dimension
        identity: Identity matrix (typically 2x2)
        
    Returns:
        The matrix representation of the specified blade
        
    Example:
        >>> I2 = np.array([[1, 0], [0, 1]])
        >>> P2 = np.array([[0, 1], [1, 0]])
        >>> # Create representation for e₁ in a 2D algebra
        >>> matrix_rep_blade(1, 2, [P2, P2], I2)
    """
    factors = []
    for j in range(dimension):
        # If the j-th bit of blade_index is 1, use the signature matrix; otherwise, use identity
        if (blade_index >> j) & 1:
            factors.append(signature_matrices[j])
        else:
            factors.append(identity)
            
    # Use the Kronecker product to combine all factors
    # The reduce function applies operator.matmul to all elements in factors
    if not factors:
        return identity
    return reduce(kronecker_product, factors)

def blade_multiplication_table(p: int = 0, q: int = 0, r: int = 0) -> np.ndarray:
    """
    Generate a multiplication table for basis blades using matrix representations.
    
    This function builds a table showing the result of multiplying any two basis blades.
    
    Args:
        p: Number of positive dimensions (default: 0)
        q: Number of negative dimensions (default: 0)
        r: Number of null dimensions (default: 0)
        
    Returns:
        A 2D numpy array where element [i,j] represents the blade index
        resulting from multiplying blade i with blade j
        
    Example:
        >>> table = blade_multiplication_table(p=2)
        >>> # Show the result of e₁ * e₂
        >>> table[1,2]
        3  # This is the index of e₁₂
    """
    d = p + q + r
    blade_count = 2**d
    
    # Generate matrix representations
    matrices = matrix_rep(p, q, r)
    
    # Initialize table
    table = np.zeros((blade_count, blade_count), dtype=np.int32)
    
    # Compute all products
    for i in range(blade_count):
        for j in range(blade_count):
            # Multiply matrices and determine the resulting blade
            product = matrices[i] @ matrices[j]
            
            # Find which basis blade this corresponds to
            # For simple cases, this can be done by comparing with all basis matrices
            for k in range(blade_count):
                # Check if product is proportional to this basis matrix
                # We need to check for both positive and negative coefficients
                if np.allclose(product, matrices[k]):
                    table[i, j] = k
                    break
                elif np.allclose(product, -matrices[k]):
                    # Store negative indices to indicate negation
                    table[i, j] = -k
                    break
            
    return table

def matrix_to_multivector(matrix: np.ndarray, matrices: List[np.ndarray]) -> Dict[int, float]:
    """
    Convert a matrix to its multivector components.
    
    This function decomposes a matrix into a linear combination of basis matrices.
    
    Args:
        matrix: The matrix to decompose
        matrices: List of basis matrices from matrix_rep()
        
    Returns:
        Dict mapping blade indices to coefficients
        
    Example:
        >>> matrices = matrix_rep(p=2)
        >>> # Create a matrix representing e₁ + 2e₂
        >>> m = matrices[1] + 2*matrices[2]
        >>> matrix_to_multivector(m, matrices)
        {1: 1.0, 2: 2.0}
    """
    blade_count = len(matrices)
    result = {}
    
    for i in range(blade_count):
        # Calculate projection of matrix onto each basis element
        # Using Frobenius inner product: <A,B> = tr(A† · B)
        basis = matrices[i]
        
        # Calculate coefficient using the inner product
        # For orthogonal basis elements, this gives the direct component
        coef = np.trace(basis.T.conjugate() @ matrix) / np.trace(basis.T.conjugate() @ basis)
        
        # Only keep non-zero coefficients with some numerical tolerance
        if abs(coef) > 1e-10:
            result[i] = float(coef.real)  # Convert to float and take real part
            
    return result