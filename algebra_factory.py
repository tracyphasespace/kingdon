"""
Algebra Factory Module for Geometric Algebra

This module provides factory methods and convenience functions 
for creating common Geometric Algebras.
"""

from __future__ import annotations

from collections import Counter
from typing import (
    Any, 
    Union, 
    Optional, 
    List, 
    Tuple, 
    Callable, 
    Dict
)

import numpy as np

from kingdon.algebra_base import AlgebraBase

# Type checking import
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from kingdon.multivector import MultiVector

class Algebra(AlgebraBase):
    """
    Full Algebra class with additional factory and convenience methods.
    """

    def multivector(self, *args: Any, **kwargs: Any) -> 'MultiVector':
        """
        Create a new multivector in this algebra.
        """
        print(f"algebra.multivector: self={type(self).__name__}, args={args}, kwargs={kwargs}")
        
        # Import here to avoid circular imports
        from kingdon.multivector import MultiVector
        
        # Directly create a MultiVector instance instead of going through MultiVectorCreation
        mv = object.__new__(MultiVector)
        mv.algebra = self
        
        # Initialize the multivector with the provided arguments
        MultiVector.__init__(mv, self, **kwargs)
        
        # Create keys and values based on the arguments
        if 'values' in kwargs and 'grades' in kwargs:
            # Get the dimension
            d = self.d
            grades = kwargs['grades']
            values = kwargs['values']
            
            # Convert to a list if it's a single value
            if not isinstance(values, (list, tuple, np.ndarray)):
                values = [values]
            
            # Determine the keys based on the grades
            if len(grades) == 1:
                grade = grades[0]
                if grade == 1:  # Vector case
                    keys = tuple(1 << i for i in range(d))
                elif grade == 2 and d >= 3 and len(values) == 3:  # Bivector case
                    keys = (3, 5, 6)
                elif grade == 0:  # Scalar case
                    keys = (0,)
                else:
                    # General case: all basis blades of the specified grade
                    keys = tuple(i for i in range(2**d) if bin(i).count('1') == grade)
            else:
                # Multiple grades: collect all basis blades of the specified grades
                keys = tuple(i for i in range(2**d) if bin(i).count('1') in grades)
            
            # Set the values and keys
            mv._keys = keys
            mv._values = values
            mv.shape = values.shape if isinstance(values, np.ndarray) and values.ndim > 1 else ()
        
        return mv

    def scalar(self, value: Union[int, float, List, Tuple, np.ndarray] = 1, 
               name: Optional[str] = None) -> 'MultiVector':
        """
        Create a scalar multivector (grade 0).
        """
        # Import here to avoid circular imports
        from kingdon.multivector import MultiVector
        
        if isinstance(value, (list, tuple, np.ndarray)):
            if len(value) != 1:
                raise ValueError(f"Expected a single scalar value, got {len(value)} values")
            value = value[0]
        
        # Create a new MultiVector with scalar grade
        mv = object.__new__(MultiVector)
        mv.algebra = self
        mv._keys = (0,)
        mv._values = [value]
        mv.shape = ()
        
        return mv

    def vector(self, values: Union[List[Any], np.ndarray, Dict[int, Any]], 
               name: Optional[str] = None) -> 'MultiVector':
        """
        Create a vector multivector (grade 1).
        
        Args:
            values: Coefficients for the vector components
            name: Optional name for symbolic representation
            
        Returns:
            A new vector multivector
        """
        # Import here to avoid circular imports
        from kingdon.multivector import MultiVector
        
        # Create a new MultiVector with vector grade
        mv = object.__new__(MultiVector)
        mv.algebra = self
        
        # Dimension determines the number of vector components
        d = self.d
        
        # Generate keys for grade-1 basis blades
        keys = tuple(1 << i for i in range(d))
        
        # Set values based on input type
        if isinstance(values, dict):
            # Dictionary input maps keys to values
            val_list = [0] * d
            for k, v in values.items():
                if 0 <= k < d:
                    val_list[k] = v
            mv._values = val_list
        else:
            # List or array input directly specifies values
            if not isinstance(values, (list, tuple, np.ndarray)):
                raise ValueError(f"Expected a list of vector components, got {type(values)}")
            if len(values) != d:
                raise ValueError(f"Expected {d} vector components, got {len(values)}")
            mv._values = list(values)
        
        mv._keys = keys
        mv.shape = ()
        
        return mv

    def bivector(self, values: Union[List[Any], np.ndarray, Dict[int, Any]], 
                 name: Optional[str] = None) -> 'MultiVector':
        """
        Create a bivector multivector (grade 2).
        
        Args:
            values: Coefficients for the bivector components
            name: Optional name for symbolic representation
            
        Returns:
            A new bivector multivector
        """
        # This method follows the same pattern as vector()
        from kingdon.multivector import MultiVector
        
        mv = object.__new__(MultiVector)
        mv.algebra = self
        
        d = self.d
        
        # Generate keys for grade-2 basis blades
        if d >= 3 and isinstance(values, (list, tuple, np.ndarray)) and len(values) == 3:
            # Special case for 3D: (e12, e13, e23) -> (3, 5, 6)
            keys = (3, 5, 6)
        else:
            # General case: all grade-2 basis blades
            keys = tuple(i for i in range(2**d) if bin(i).count('1') == 2)
        
        # Set values
        if isinstance(values, dict):
            val_list = [0] * len(keys)
            for i, k in enumerate(keys):
                val_list[i] = values.get(k, 0)
            mv._values = val_list
        else:
            if len(values) != len(keys):
                raise ValueError(f"Expected {len(keys)} bivector components, got {len(values)}")
            mv._values = list(values)
        
        mv._keys = keys
        mv.shape = ()
        
        return mv

    def grade(self, mv: 'MultiVector', *grades: int) -> 'MultiVector':
        """
        Extract specified grades from a multivector.
        
        Args:
            mv: Input multivector
            *grades: Grades to extract
            
        Returns:
            A new multivector containing only the specified grades
        """
        return mv.grade(*grades)

    @classmethod
    def from_signature(cls, signature: List[int]) -> 'Algebra':
        """
        Create an algebra from a custom signature.
        
        Args:
            signature: List of 1, -1, and 0 values defining the signature
            
        Returns:
            A new Algebra instance
        """
        counts = Counter(signature)
        p, q, r = counts[1], counts[-1], counts[0]
        return cls(p=p, q=q, r=r, signature=signature)

    @classmethod
    def from_name(cls, name: str) -> 'Algebra':
        """
        Create a common algebra by name.
        
        Args:
            name: Name of the algebra to create
            
        Returns:
            A new Algebra instance
            
        Raises:
            ValueError: If the algebra name is not recognized
        """
        name = name.upper()
        
        algebras = {
            "CGA2D": (3, 1, 0),   # 2D Conformal Geometric Algebra
            "CGA3D": (4, 1, 0),   # 3D Conformal Geometric Algebra
            "PGA2D": (2, 0, 1),   # 2D Projective Geometric Algebra
            "PGA3D": (3, 0, 1),   # 3D Projective Geometric Algebra
            "VGA2D": (2, 0, 0),   # 2D Vector Geometric Algebra
            "VGA3D": (3, 0, 0),   # 3D Vector Geometric Algebra
            "STA": (1, 3, 0),     # Spacetime Algebra
            "APS": (0, 3, 0)      # Algebra of Physical Space
        }
        
        if name not in algebras:
            raise ValueError(f"Unknown algebra name: {name}")
        
        p, q, r = algebras[name]
        return cls(p=p, q=q, r=r)

# Convenience factory functions
def PGA2D() -> Algebra:
    """Create a 2D Projective Geometric Algebra."""
    return Algebra(p=2, q=0, r=1)

def PGA3D() -> Algebra:
    """Create a 3D Projective Geometric Algebra."""
    return Algebra(p=3, q=0, r=1)

def CGA2D() -> Algebra:
    """Create a 2D Conformal Geometric Algebra."""
    return Algebra(p=3, q=1, r=0)

def CGA3D() -> Algebra:
    """Create a 3D Conformal Geometric Algebra."""
    return Algebra(p=4, q=1, r=0)

def VGA3D() -> Algebra:
    """Create a 3D Vector Geometric Algebra."""
    return Algebra(p=3, q=0, r=0)

def STA() -> Algebra:
    """Create a Spacetime Algebra."""
    return Algebra(p=1, q=3, r=0)

def APS() -> Algebra:
    """Create an Algebra of Physical Space."""
    return Algebra(p=0, q=3, r=0)