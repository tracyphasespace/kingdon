"""
Kingdon: A Geometric Algebra Library for Python
===============================================

Kingdon is a Python library for Geometric Algebra (GA) computations. It provides
an efficient implementation of Geometric Algebra suitable for physics, graphics,
and robotics.

Basic Usage:
-----------
>>> from kingdon import Algebra
>>> alg = Algebra(p=3, q=0, r=0)
>>> from kingdon.ga_factory import create_vector
>>> v = create_vector(alg, [1, 0, 0])
>>> print(v)
"""

__author__ = "Martin Roelfs, Modified by Tracy and David"
__email__ = 'martinroelfs@yahoo.com'
__version__ = '1.3.1'

from sympy import symbols

# Core components
from kingdon.algebra import Algebra
from kingdon.multivector import MultiVector
from .operator_dict import AlgebraError # <--- ADDED THIS LINE

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

# Visualization (optional)
try:
    from kingdon.graph import GraphWidget
except ImportError:
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

# Optional: Define __all__ to control `from kingdon import *`
__all__ = [
    'Algebra',
    'MultiVector',
    'AlgebraError', # Add AlgebraError here too
    'create_vector',
    'create_bivector',
    'create_scalar',
    'create_trivector',
    'create_pseudoscalar',
    'create_blade',
    'matrix_rep',
    'PGA2D', 'PGA3D', 'CGA2D', 'CGA3D', 'VGA3D', 'STA',
    'version',
    'symbols' # Re-export symbols if desired
]

# Conditionally add GraphWidget to __all__ if imported
if 'GraphWidget' in locals():
    __all__.append('GraphWidget')