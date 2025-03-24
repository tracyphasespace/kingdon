from __future__ import annotations

"""
MultiVector Creation Module

Handles the creation and initialization of MultiVector objects with 
various input types and configurations.
"""

# Standard library and typing imports
import operator
import re
from typing import (
    Any, 
    Callable, 
    Dict, 
    Optional, 
    Tuple, 
    Union, 
    List, 
    TypeVar,
    Type
)

# Third-party imports
import numpy as np
import sympy
from sympy import Expr

# Local imports
from kingdon.multivector_base import MultiVectorBase
try:
    from kingdon.polynomial import RationalPolynomial
except ImportError:
    class RationalPolynomial:
        pass

# Type variable for scalar values
T = TypeVar('T', int, float, complex, Expr, RationalPolynomial)

class MultiVectorCreation(MultiVectorBase):
    """
    Extension of MultiVectorBase focused on object creation and initialization.
    """
    
    @classmethod
    def fromkeysvalues(cls, algebra: Any, keys: Tuple[int, ...], values: Union[List[Any], np.ndarray]) -> 'MultiVectorCreation':
        """
        Construct a multivector from provided keys and values.
        """
        print(f"MultiVectorCreation.fromkeysvalues: algebra={type(algebra).__name__}, keys={keys}, values={values}")
        obj = object.__new__(cls)
        obj.algebra = algebra
        obj._keys = keys
        obj._values = values if isinstance(values, (np.ndarray, list)) else [values]
        obj.shape = values.shape if isinstance(values, np.ndarray) and values.ndim > 1 else ()
        return obj
    
    @classmethod
    def __new__(cls, algebra: Any, *args, **kwargs) -> 'MultiVectorCreation':
        """
        Create a new MultiVector with flexible initialization methods.
        """
        print(f"MultiVectorCreation.__new__: cls={cls.__name__}, algebra={type(algebra).__name__}, args={args}, kwargs={kwargs}")
        
        # Extract kwargs with defaults
        values = kwargs.get('values', None)
        keys = kwargs.get('keys', None)
        name = kwargs.get('name', None)
        grades = kwargs.get('grades', None)
        symbolcls = kwargs.get('symbolcls', sympy.Symbol)

        # Handle positional values argument
        if args and values is None:
            values = args[0]
            args = args[1:]

        # Ensure we access the d attribute correctly - crucial for fixing this bug
        # in this context, algebra appears to be a class, not an instance
        try:
            # When a proper instance is passed
            d = algebra.d
        except AttributeError:
            # Access d directly from the first argument - which contains the actual algebra instance
            if args and hasattr(args[0], 'd'):
                d = args[0].d
            else:
                # Default dimension if all else fails
                d = 3
                print(f"Warning: Could not determine algebra dimension. Using default d={d}")

        # Handle named blade initialization from kwargs
        blade_kwargs = {k: v for k, v in kwargs.items() if k not in ('keys', 'values', 'name', 'grades', 'symbolcls')}
        if blade_kwargs and not values and not keys:
            processed_kwargs = {}
            # Use getattr to safely access the canon2bin attribute
            canon2bin = getattr(algebra, 'canon2bin', {})
            for blade_name, value in blade_kwargs.items():
                blade_key = canon2bin.get(blade_name, blade_name)
                processed_kwargs[blade_key] = value
            keys = list(processed_kwargs.keys())
            values = list(processed_kwargs.values())

        # Special case for scalar creation
        if grades == (0,) and values is not None:
            if not isinstance(values, (list, tuple, np.ndarray)):
                values = [values]
            elif isinstance(values, (list, tuple, np.ndarray)) and len(values) == 1:
                values = list(values)
            return cls.fromkeysvalues(algebra, (0,), values)

        # Determine grades if not provided
        if grades is None:
            if keys:
                grades = tuple(sorted(set(bin(k).count('1') for k in keys if isinstance(k, int))))
            elif name:
                grades = tuple(range(d + 1))
            else:
                grades = (0,)

        # Handle specific grade cases
        if values is not None and keys is None:
            if len(grades) == 1:
                grade = grades[0]
                if grade == 1 and d == len(values):  # Vector case
                    keys = tuple(1 << i for i in range(d))
                    values = [v for v in values]  # Copy values
                elif grade == 2 and len(values) == 3 and d >= 3:  # Bivector case (minimum 3D)
                    keys = (3, 5, 6)
                    values = [v for v in values]  # Copy values
                elif grade == 0:  # Scalar case
                    keys = (0,)
                    values = [values] if not isinstance(values, (list, tuple, np.ndarray)) else list(values)
                else:
                    keys = tuple(i for i in range(2**d) if bin(i).count('1') == grade)
                    if not keys:
                        raise ValueError(f"No basis blades exist for grade {grade} in {d}-dimensional algebra")
                    if len(values) == 1:
                        values = [values[0]] * len(keys)  # Broadcast scalar
                    elif len(values) != len(keys):
                        raise ValueError(f"Expected {len(keys)} values for grade {grade} in {d}-dimensional algebra, got {len(values)}")

        # Default to empty multivector if no values
        if values is None:
            if keys is None:
                keys = tuple(i for i in range(2**d) if bin(i).count('1') in grades)
            values = [0] * len(keys)

        # Create symbolic values if name provided
        if name and all(v == 0 for v in values):
            values = [symbolcls(f"{name}_{k}") for k in keys]

        # Final validation
        if not (isinstance(values, np.ndarray) and values.ndim > 1) and len(keys) != len(values):
            raise ValueError(f"Length of `keys` ({len(keys)}) and `values` ({len(values)}) must match.")

        return cls.fromkeysvalues(algebra, tuple(keys), values)

# Export the MultiVectorCreation class
__all__ = ['MultiVectorCreation']