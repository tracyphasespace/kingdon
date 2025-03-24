"""
Kingdon: A Geometric Algebra Library for Python
===============================================

Kingdon is a Python library for Geometric Algebra (GA) computations. It provides
a flexible and high-performance implementation of Geometric Algebra suitable for
use in physics, computer graphics, robotics, and other fields.

Key Components:
--------------
* Algebra: Create and work with different Geometric Algebras
* MultiVector: Manipulate elements of Geometric Algebras
* Factory functions: Easily create common GA objects
* Matrix representations: Convert GA operations to matrix form
* Visualization tools: Visualize geometric objects

Basic Usage:
-----------
>>> from kingdon import Algebra
>>> from kingdon.ga_factory import create_vector
>>> 
>>> # Create a 3D Euclidean Geometric Algebra
>>> alg = Algebra(p=3, q=0, r=0)
>>> 
>>> # Create two vectors
>>> v1 = create_vector(alg, [1, 0, 0])
>>> v2 = create_vector(alg, [0, 1, 0])
>>> 
>>> # Compute the outer product (bivector)
>>> B = v1 ^ v2
>>> print(B)
1𝐞₁₂
"""

__author__ = """Martin Roelfs, Modified by Tracy and David """
__email__ = 'martinroelfs@yahoo.com'
__version__ = '1.3.0'

# Core components
from kingdon.algebra import Algebra
from kingdon.multivector import MultiVector

# Factory functions for easy object creation
from kingdon.ga_factory import (
    create_vector,
    create_bivector,
    create_scalar,
    create_trivector,
    create_pseudoscalar,
    create_blade
)

# Matrix representation functions
from kingdon.matrixreps import matrix_rep

# Visualization
try:
    from kingdon.graph import GraphWidget
except ImportError:
    # Optional visualization component
    pass

# Convenient abbreviations for common algebras
def PGA2D() -> Algebra:
    """Create a 2D Projective Geometric Algebra (p=2, q=0, r=1)."""
    return Algebra(p=2, q=0, r=1)

def PGA3D() -> Algebra:
    """Create a 3D Projective Geometric Algebra (p=3, q=0, r=1)."""
    return Algebra(p=3, q=0, r=1)

def CGA2D() -> Algebra:
    """Create a 2D Conformal Geometric Algebra (p=3, q=1, r=0)."""
    return Algebra(p=3, q=1, r=0)

def CGA3D() -> Algebra:
    """Create a 3D Conformal Geometric Algebra (p=4, q=1, r=0)."""
    return Algebra(p=4, q=1, r=0)

def VGA3D() -> Algebra:
    """Create a 3D Vector Geometric Algebra (p=3, q=0, r=0)."""
    return Algebra(p=3, q=0, r=0)

def STA() -> Algebra:
    """Create a Spacetime Algebra (p=1, q=3, r=0)."""
    return Algebra(p=1, q=3, r=0)

# Version information
def version() -> str:
    """Return the version of the Kingdon library."""
    return __version__