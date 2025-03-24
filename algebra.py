"""
Algebra Module for Geometric Algebra

This module provides a unified import point for Algebra-related classes and functions.
"""

from __future__ import annotations

# Import the base Algebra components
from kingdon.algebra_base import AlgebraBase

# Import the factory Algebra and convenience functions
from kingdon.algebra_factory import (
    Algebra, 
    PGA2D, 
    PGA3D, 
    CGA2D, 
    CGA3D, 
    VGA3D, 
    STA, 
    APS
)

# Re-export everything for easy importing
__all__ = [
    'AlgebraBase',
    'Algebra',
    'PGA2D',
    'PGA3D',
    'CGA2D',
    'CGA3D',
    'VGA3D',
    'STA',
    'APS'
]