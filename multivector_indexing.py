"""
MultiVector Indexing and Assignment Module

Provides advanced indexing and assignment capabilities for MultiVector objects.
"""

from __future__ import annotations
import copy
import re
from typing import (
    Any, Callable, Dict, Optional, 
    Tuple, Union, TypeVar, List, cast
)

import numpy as np
import sympy

# Conditional import to avoid circular dependency
try:
    from kingdon.multivector_creation import MultiVectorCreation
except ImportError:
    # Fallback if direct import fails
    MultiVectorCreation = object

class MultiVectorIndexing(MultiVectorCreation):
    """
    Extension of MultiVectorCreation with advanced indexing and assignment methods.
    
    This class provides methods for indexing into multivectors and assigning values,
    with special handling for NumPy arrays and other array-like structures.
    """
    
    def __getitem__(self, item: Union[int, Tuple, slice]) -> MultiVectorIndexing:
        """
        Allow flexible indexing into the multivector's coefficients.
        
        Supports:
        - Single integer index
        - Tuple indices
        - Slicing
        - Extracting array-valued multivectors
        """
        # Make sure we have the required attributes
        if not hasattr(self, '_keys'):
            self._keys = getattr(self, 'keys', lambda: (0,))()
            
        if not hasattr(self, '_values'):
            self._values = getattr(self, 'values', lambda: [0])()
            
        if not hasattr(self, 'shape'):
            self.shape = ()
    
        if not isinstance(item, tuple):
            item = (item,)
    
        vals = self.values()
    
        # Handle array-like values
        if isinstance(vals, (list, tuple)):
            # For single integer index, return that component
            if len(item) == 1 and isinstance(item[0], int):
                return self.fromkeysvalues(
                    self.algebra, 
                    (self._keys[item[0]],), 
                    [vals[item[0]]]
                )
        
            # For other indexing, try to apply to each component
            new_vals = []
            for v in vals:
                if hasattr(v, '__getitem__'):
                    try:
                        new_vals.append(v[item])
                    except (TypeError, IndexError):
                        new_vals.append(v)
                else:
                    new_vals.append(v)
        elif isinstance(vals, np.ndarray):
            if vals.ndim > 1:
                # For multidimensional arrays, apply the index
                new_vals = vals[item]
                # Important: Preserve shape information for sliced arrays
                if isinstance(new_vals, np.ndarray) and new_vals.ndim == 1:
                    # If we sliced to a 1D array, make sure to return a MultiVector 
                    # with the proper shape
                    mv = self.fromkeysvalues(
                        self.algebra, 
                        self.keys(), 
                        new_vals
                    )
                    mv.shape = (len(new_vals),)  # Explicitly set the shape
                    return mv
            else:
                # For 1D arrays, handle integer indexing
                if len(item) == 1 and isinstance(item[0], int):
                    return self.fromkeysvalues(
                        self.algebra, 
                        (self._keys[item[0]],), 
                        [vals[item[0]]]
                    )
                new_vals = vals[item]
        else:
            # For scalar values, just return them
            new_vals = [vals]
            
        return self.fromkeysvalues(
            self.algebra, 
            self.keys(), 
            new_vals
        )

    def __setitem__(self, key: Union[int, Tuple, slice], value: Any) -> None:
        """
        Handle assignment to MultiVector objects with improved NumPy interoperability.
        
        Supports:
        - Assigning values to specific components
        - Assigning values to array-valued multivectors
        - Assigning values from NumPy arrays
        """
        # Make sure we have the required attributes
        if not hasattr(self, '_keys'):
            self._keys = getattr(self, 'keys', lambda: (0,))()
            
        if not hasattr(self, '_values'):
            self._values = getattr(self, 'values', lambda: [0])()
    
        # Validate key type
        if not isinstance(key, tuple):
            key = (key,)
        
        # Handle MultiVector assignment
        if hasattr(value, 'values') and callable(value.values):
            # For a single MultiVector, assign its values
            if len(value.values()) == 1:
                value = value.values()[0]
            else:
                value = value.values()
        
        # Handle NumPy array assignment
        if isinstance(value, np.ndarray):
            # Ensure self._values is a NumPy array for consistent handling
            if not isinstance(self._values, np.ndarray):
                self._values = np.array(self._values)
            
            # For object dtype arrays of MultiVectors, handle special case
            if value.dtype == object:
                # Replace individual elements with their values or the whole MultiVector
                for i, v in np.ndenumerate(value):
                    if hasattr(v, 'values') and callable(v.values):
                        # If single-valued, use its first value
                        if len(v.values()) == 1:
                            self._values[key + i] = v.values()[0]
                        else:
                            self._values[key + i] = v.values()
            else:
                # For regular numeric arrays, assign normally
                self._values[key] = value
        
        # Handle scalar assignments and list assignments
        elif isinstance(self._values, np.ndarray):
            # For array-valued multivectors, assign using NumPy
            self._values[key] = value
        else:
            # For list-based multivectors, fall back to standard list assignment
            if isinstance(self._values, list):
                self._values[key[0] if len(key) == 1 else key] = value
            else:
                # Fallback for other unexpected cases
                self._values[key[0] if len(key) == 1 else key] = value

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        """
        Support for NumPy's universal functions.
        
        This method enables multivectors to work with NumPy's mathematical functions.
        """
        # Extract the multivector values and apply the ufunc
        if method == '__call__':
            # Extract values from multivector inputs
            args = []
            for input_arg in inputs:
                if isinstance(input_arg, MultiVectorIndexing):
                    args.append(input_arg.values())
                else:
                    args.append(input_arg)
                    
            # Apply the ufunc to the values
            result = ufunc(*args, **kwargs)
            
            # Wrap the result in a new multivector
            return self.__class__.fromkeysvalues(self.algebra, self._keys, result)
            
        return NotImplemented