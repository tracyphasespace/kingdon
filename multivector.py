"""
MultiVector Module for Geometric Algebra

This module consolidates the MultiVector class from its component modules.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

# Import all required components
from kingdon.multivector_base import MultiVectorBase
from kingdon.multivector_creation import MultiVectorCreation
from kingdon.multivector_indexing import MultiVectorIndexing
from kingdon.multivector_operations import MultiVectorOperations
from kingdon.multivector_transforms import MultiVectorTransforms
from kingdon.multivector_core import MultiVector as CoreMultiVector

# For type checking only
if TYPE_CHECKING:
    try:
        from kingdon.algebra import Algebra
    except ImportError:
        # During initial import
        Algebra = object

# Define the main MultiVector class that inherits from all components
class MultiVector(MultiVectorTransforms, MultiVectorOperations, CoreMultiVector, MultiVectorIndexing, MultiVectorCreation, MultiVectorBase):
    """
    Consolidated MultiVector class that includes all functionality from component classes.
    
    This class provides a complete implementation of multivectors for geometric algebra,
    including creation, indexing, operations, and transformations.
    """
    pass

# Export the MultiVector class
__all__ = ['MultiVector']