"""
Utility Functions for Geometric Algebra
======================================

This module provides utility functions for working with MultiVector objects 
and implementing common geometric algebra operations.
"""

import re
import math
import numpy as np
from typing import Any, List, Tuple, Dict, Callable, Optional, Union, Generator

def map_multivector(mv, func):
    """
    Apply a function to each component of a multivector.
    
    Args:
        mv: The input multivector
        func: Function to apply to each component value
        
    Returns:
        A new multivector with transformed values
    """
    # Apply the function to each value
    new_values = [func(v) for v in mv._values]
    
    # Create a new multivector with the transformed values
    return mv.__class__.fromkeysvalues(mv.algebra, mv._keys, new_values)

def filter_multivector(mv, pred):
    """
    Filter components of a multivector based on a predicate function.
    
    Args:
        mv: The input multivector
        pred: Predicate function that takes a value and returns True/False
        
    Returns:
        A new multivector with only the components for which pred returns True
    """
    # Filter components
    filtered_items = [(k, v) for k, v in zip(mv._keys, mv._values) if pred(v)]
    
    if not filtered_items:
        return mv.__class__.fromkeysvalues(mv.algebra, (), [])
        
    keys, values = zip(*filtered_items)
    return mv.__class__.fromkeysvalues(mv.algebra, keys, list(values))

def exp(mv):
    """
    Compute the exponential of a multivector.
    
    Args:
        mv: Input multivector
        
    Returns:
        A new multivector representing exp(mv)
    """
    # Special case for zero
    if not mv._keys or all(v == 0 for v in mv._values):
        return mv.algebra.scalar(1)
        
    # For vectors and bivectors, use specialized formulas
    if len(mv.grades) == 1 and mv.grades[0] in (1, 2):
        # Square the multivector
        mv2 = mv * mv
        
        # Check if it's nilpotent (x² = 0)
        if mv2.grades == (0,) and abs(mv2.e) < 1e-10:
            return mv.algebra.scalar(1) + mv
            
        # Check if it has a negative square (typical for bivectors)
        if mv2.grades == (0,) and mv2.e < 0:
            # Use formula exp(x) = cos(|x|) + x/|x| * sin(|x|)
            alpha = math.sqrt(-mv2.e)
            return mv.algebra.scalar(math.cos(alpha)) + (mv / alpha) * math.sin(alpha)
            
        # If it has a positive square, use hyperbolic functions
        if mv2.grades == (0,) and mv2.e > 0:
            alpha = math.sqrt(mv2.e)
            return mv.algebra.scalar(math.cosh(alpha)) + (mv / alpha) * math.sinh(alpha)
            
    # For general multivectors, use Taylor series
    result = mv.algebra.scalar(1)
    term = mv.algebra.scalar(1)
    factorial = 1
    
    # Compute the series up to 15 terms
    for i in range(1, 16):
        factorial *= i
        term = term * mv
        result += term / factorial
        
        # Stop if the term becomes negligible
        if all(abs(v) < 1e-10 for v in term._values):
            break
            
    return result

def outersin(mv):
    """
    Compute the outer sine of a multivector.
    
    Args:
        mv: Input multivector
        
    Returns:
        A new multivector representing sin^(mv)
    """
    # Special case for zero
    if not mv._keys or all(v == 0 for v in mv._values):
        return mv.algebra.scalar(0)
        
    # For grade-2 bivector, use specialized formula
    if mv.grades == (2,):
        # Square the bivector
        mv2 = mv * mv
        
        # Check if it has a negative square (typical for bivectors)
        if mv2.grades == (0,) and mv2.e < 0:
            # Use formula sin^(B) = B/|B| * sin(|B|)
            alpha = math.sqrt(-mv2.e)
            return (mv / alpha) * math.sin(alpha)
            
    # For general multivectors, use Taylor series with only odd terms
    result = mv.algebra.scalar(0)
    term = mv
    
    # Add the first term
    result += term
    
    # Compute the series up to 8 terms (only odd powers)
    for i in range(3, 16, 2):
        factorial = math.factorial(i)
        term = term * mv * mv / ((i-1) * i)
        
        # Alternate signs
        if (i-1) // 2 % 2 == 1:
            result -= term
        else:
            result += term
        
        # Stop if the term becomes negligible
        if all(abs(v) < 1e-10 for v in term._values):
            break
            
    return result

def outercos(mv):
    """
    Compute the outer cosine of a multivector.
    
    Args:
        mv: Input multivector
        
    Returns:
        A new multivector representing cos^(mv)
    """
    # Special case for zero
    if not mv._keys or all(v == 0 for v in mv._values):
        return mv.algebra.scalar(1)
        
    # For grade-2 bivector, use specialized formula
    if mv.grades == (2,):
        # Square the bivector
        mv2 = mv * mv
        
        # Check if it has a negative square (typical for bivectors)
        if mv2.grades == (0,) and mv2.e < 0:
            # Use formula cos^(B) = cos(|B|)
            alpha = math.sqrt(-mv2.e)
            return mv.algebra.scalar(math.cos(alpha))
            
    # For general multivectors, use Taylor series with only even terms
    result = mv.algebra.scalar(1)
    term = mv.algebra.scalar(1)
    
    # Compute the series up to 8 terms (only even powers)
    for i in range(2, 17, 2):
        factorial = math.factorial(i)
        term = term * mv * mv / ((i-1) * i)
        
        # Alternate signs
        if (i // 2) % 2 == 1:
            result -= term
        else:
            result += term
        
        # Stop if the term becomes negligible
        if all(abs(v) < 1e-10 for v in term._values):
            break
            
    return result

def outertan(mv):
    """
    Compute the outer tangent of a multivector.
    
    Args:
        mv: Input multivector
        
    Returns:
        A new multivector representing tan^(mv)
    """
    # Compute using the ratio of outer sine to outer cosine
    cos_mv = outercos(mv)
    sin_mv = outersin(mv)
    
    # Check if cosine is near zero
    cosine_norm = cos_mv.norm()
    if cosine_norm < 1e-10:
        raise ValueError("Outer tangent undefined: cosine is zero")
        
    return sin_mv / cos_mv

def create_reciprocal_frame(vectors):
    """
    Compute the reciprocal frame for a set of vectors.
    
    Args:
        vectors: List of vectors that form a frame
        
    Returns:
        List of reciprocal frame vectors
    """
    if not vectors:
        return []
        
    n = len(vectors)
    algebra = vectors[0].algebra
    
    # Compute the pseudoscalar for the frame
    pseudoscalar = vectors[0]
    for i in range(1, n):
        pseudoscalar = pseudoscalar ^ vectors[i]
        
    # Compute the reciprocal frame vectors
    reciprocal = []
    for i in range(n):
        # Remove the i-th vector from the wedge product
        blade = algebra.scalar(1)
        for j in range(n):
            if j != i:
                blade = blade ^ vectors[j]
                
        # Compute the reciprocal frame vector
        r_vector = (-1)**(i) * blade.dual() / pseudoscalar
        reciprocal.append(r_vector)
        
    return reciprocal

def rotate_vector(rotor, vector):
    """
    Rotate a vector using a rotor.
    
    Args:
        rotor: The rotor (must be normalized)
        vector: The vector to rotate
        
    Returns:
        Rotated vector
    """
    # Apply the sandwich product: R * v * ~R
    return rotor * vector * ~rotor

def reflect_vector(normal, vector):
    """
    Reflect a vector in a hyperplane with the given normal vector.
    
    Args:
        normal: Normal vector to the reflection plane (should be normalized)
        vector: The vector to reflect
        
    Returns:
        Reflected vector
    """
    # Apply the formula: -n * v * n
    return -normal * vector * normal

def project_vector(vector, onto):
    """
    Project a vector onto another vector.
    
    Args:
        vector: The vector to project
        onto: The vector to project onto
        
    Returns:
        Projected vector
    """
    # Compute the inner product
    inner = vector | onto
    
    # Compute the inverse of the "onto" vector
    onto_inv = onto / onto.norm()**2
    
    # Project using the formula: (v·a) * a⁻¹
    return inner * onto_inv

def matrix_to_multivector(matrix, matrices, algebra):
    """
    Convert a matrix to its multivector representation.

    Args:
        matrix: The matrix to convert (2D NumPy array)
        matrices: List of basis matrices from matrix_rep()
        algebra: The Geometric Algebra instance

    Returns:
        MultiVector representing the matrix

    Raises:
        ValueError: if shapes mismatch or inputs are invalid
        TypeError: if algebra has no dimension attribute
    """
    import warnings
    from .multivector import MultiVector

    # Validate inputs
    if not hasattr(algebra, 'd'):
        raise TypeError(f"Invalid algebra: missing 'd' attribute ({type(algebra)})")
    if not isinstance(matrices, (list, tuple)) or len(matrices) == 0:
        raise ValueError("`matrices` must be a non-empty list of basis matrices")
    if not hasattr(matrix, 'shape') or matrix.shape != matrices[0].shape:
        raise ValueError(f"Shape mismatch: matrix {getattr(matrix, 'shape', None)} vs basis {matrices[0].shape}")

    coefs: Dict[int, float] = {}
    for i, basis in enumerate(matrices):
        denom = np.sum(basis * basis.conj().T)
        if abs(denom) < 1e-16:
            warnings.warn(f"Skipping basis index {i}: zero Frobenius norm")
            continue

        coef = np.sum(matrix * basis.conj().T) / denom
        if abs(coef) > 1e-10:
            key = 0 if i == 0 else (1 << (i - 1))
            coefs[key] = float(coef.real)

    if not coefs:
        return algebra.scalar(0)

    keys = tuple(coefs.keys())
    values = [coefs[k] for k in keys]
    return MultiVector.fromkeysvalues(algebra, keys, values)


def swap_blades(blade_str):
    """
    Compute the canonical form of a blade string and the number of swaps needed.
    
    Args:
        blade_str: String representation of a blade (e.g., 'e21')
        
    Returns:
        Tuple (canonical_form, swaps) where canonical_form is the sorted blade string
        and swaps is the number of swaps needed to sort it
    """
    # Handle special cases
    if blade_str in ("e", "e0"):
        return blade_str, 0
        
    # Extract indices
    match = re.match(r'^e([0-9a-fA-F]+)$', blade_str)
    if not match:
        raise ValueError(f"Invalid basis blade format: {blade_str}")
        
    # Get indices and count swaps using bubble sort
    indices = [int(c) for c in match.group(1)]
    sorted_indices = sorted(indices)
    
    swaps = 0
    # Count swaps using bubble sort
    for i in range(len(indices)):
        for j in range(len(indices) - 1, i, -1):
            if indices[j-1] > indices[j]:
                indices[j-1], indices[j] = indices[j], indices[j-1]
                swaps += 1
                
    # Create canonical form
    canonical = f"e{''.join(str(idx) for idx in sorted_indices)}"
    
    return canonical, swaps