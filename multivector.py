# Created by install_patches.py
"""
MultiVector Module for Geometric Algebra

This module provides the MultiVector class, which represents elements in 
Geometric Algebra (GA) and their operations.
"""
# At the top of multivector.py, replace the import section with:

from __future__ import annotations

import copy
import re
import string
import operator
from functools import reduce
from typing import (
    Any, Callable, Dict, Generator, List, Optional, Sequence,
    Set, Tuple, Union, cast, TYPE_CHECKING
)

import numpy as np
import sympy
from sympy import Expr

# Use TYPE_CHECKING for imports that might cause circular imports
if TYPE_CHECKING:
    from kingdon.polynomial import RationalPolynomial
else:
    # Create a placeholder class for runtime checks
    class RationalPolynomial:
        """Placeholder for RationalPolynomial to avoid circular imports."""
        pass

class MultiVector:
    """
    Comprehensive MultiVector class for Geometric Algebra (GA).
    
    A multivector represents an element in a Geometric Algebra, which can be
    a combination of different geometric entities (scalar, vector, bivector, etc.).
    
    Attributes:
        algebra: The Geometric Algebra instance this multivector belongs to
        _keys: Tuple of binary keys corresponding to basis blades
        _values: Coefficients for each basis blade (list or numpy array)
        shape: Shape of the multivector (empty tuple for scalar, or specific shape for array-valued)
    """
    
    def get(self, key: int, default: Any = None) -> Any:
        """
        Get the coefficient for a given basis blade key, with a default if not found.
        
        Args:
            key: Binary key of the basis blade
            default: Value to return if key is not found
            
        Returns:
            The coefficient associated with the key, or default if not present
        """
        try:
            idx = self._keys.index(key)
            return self._values[idx]
        except (ValueError, IndexError):
            return default
    
    __array_priority__ = -100

    def __array__(self, dtype=None):
        """
        Support conversion to NumPy array.
        
        Args:
            dtype: Data type for the array
                
        Returns:
            NumPy array of the values
        """
        if dtype == object:
            # For object arrays, we need to return a numpy array with self as an element
            # The key is to create a new object array that *contains* self, rather than returning self directly
            arr = np.empty(1, dtype=object)
            arr[0] = self
            return arr
            
        # For numeric arrays, extract values
        values = list(self._values)
        if not values:
            values = [0]
        
        return np.array(values, dtype=dtype)
    
    
    def __new__(cls, *args, **kwargs) -> MultiVector:
        """
        Create a new MultiVector instance.
        
        Args:
            *args: Positional arguments, with the first being the algebra
            **kwargs: Keyword arguments for initialization
            
        Returns:
            A new MultiVector instance
        """
        
        # Create a new instance
        obj = object.__new__(cls)
        
        # Return the new instance
        return obj
    
    def __init__(self, algebra: Any, **kwargs: Any) -> None:
        """
        Initialize a MultiVector with provided parameters.
        
        Args:
            algebra: The Geometric Algebra instance
            **kwargs: Initialization parameters including:
                - keys: Binary keys for basis blades
                - values: Coefficients for basis blades
                - name: Optional name for symbolic representation
                - grades: Filter to specific grades
        """
        
        # Initialize core attributes
        self.algebra = algebra
        
        # Only set _keys and _values if they haven't been set already
        if not hasattr(self, '_keys'):
            self._keys = kwargs.get('keys', (0,))
            
        if not hasattr(self, '_values'):
            if 'values' in kwargs:
                values = kwargs['values']
                self._values = values if isinstance(values, (list, np.ndarray)) else [values]
            else:
                self._values = [0] * len(self._keys)
        
        if not hasattr(self, 'shape'):
            self.shape = getattr(self._values, 'shape', ()) if isinstance(self._values, np.ndarray) and self._values.ndim > 1 else ()
    
    @classmethod
    def fromkeysvalues(cls, algebra: Any, keys: Tuple[int, ...], values: Union[List[Any], np.ndarray]) -> 'MultiVector':
        """
        Construct a multivector directly from keys and values.
        
        Args:
            algebra: The Geometric Algebra instance
            keys: Binary keys for basis blades
            values: Coefficients for the basis blades
                
        Returns:
            A new MultiVector instance
        """
        # Ensure keys is always a tuple, even if a single integer is passed
        if isinstance(keys, int):
            keys = (keys,)
        
        obj = object.__new__(cls)
        obj.algebra = algebra
        obj._keys = keys
        obj._values = values
        obj.shape = values.shape if isinstance(values, np.ndarray) and values.ndim > 1 else ()
        return obj
    
    # === Properties and Attribute Access ===
    
    @property
    def d(self) -> int:
        """
        Get the dimension of the algebra.
        
        Returns:
            Total number of dimensions in the algebra
        """
        return self.algebra.d if hasattr(self.algebra, 'd') else 0
    
    @property
    def coeffs(self) -> Dict[int, Any]:
        """
        Get a dictionary mapping basis blade keys to their coefficients.
        
        For array-valued multivectors, returns the coefficients for the first element.
        
        Returns:
            Dictionary mapping binary keys to their coefficients
        """
        if isinstance(self._values, np.ndarray) and self._values.ndim > 1:
            # For array-valued multivectors, return the first element's coefficients
            return dict(zip(self._keys, self._values.flat[:len(self._keys)]))
        return dict(zip(self._keys, self._values))
    
    @property
    def grades(self) -> Tuple[int, ...]:
        """
        Get the sorted tuple of grades present in the multivector.
        
        Returns:
            Sorted tuple of grades (0=scalar, 1=vector, 2=bivector, etc.)
        """
        if not hasattr(self, '_keys') or not self._keys:
            return tuple()
                
        # Ensure _keys is iterable
        keys = self._keys
        if isinstance(keys, int):
            keys = (keys,)  # Convert integer to tuple with one element
                
        # Filter out non-zero coefficients
        grades_dict = {}
        for k, v in zip(keys, self._values):
            grade = bin(k).count("1")
            # Only include non-zero components
            is_zero = (
                (isinstance(v, (int, float, complex)) and v == 0) or
                (isinstance(v, np.ndarray) and np.all(v == 0))
            )
            if not is_zero:
                grades_dict[grade] = True
                    
        # If no non-zero coefficients found, ensure we still return grade 0 for zero scalar
        if not grades_dict and hasattr(self, 'keys') and callable(self.keys):
            keys_list = list(self.keys())
            if keys_list and 0 in [bin(k).count("1") for k in keys_list]:
                grades_dict[0] = True
                
        return tuple(sorted(grades_dict.keys()))
    
    
    
    
    @property
    def issymbolic(self) -> bool:
        """
        Determine if the multivector contains symbolic coefficients.
        
        Returns:
            True if any coefficient is symbolic, False otherwise
        """
        # Import locally to avoid circular imports
        from kingdon.polynomial import RationalPolynomial, Polynomial
        
        # Check each value in the multivector
        for v in self.values():
            # Check for sympy Expr objects
            if isinstance(v, sympy.Expr) and not v.is_number:
                return True
                
            # Check for objects with free_symbols
            if hasattr(v, "free_symbols") and v.free_symbols:
                return True
                
            # Check for Polynomial or RationalPolynomial
            if isinstance(v, (Polynomial, RationalPolynomial)):
                return True
                
            # Check if value is a custom symbolic class from algebra
            if hasattr(self.algebra, 'codegen_symbolcls') and self.algebra.codegen_symbolcls:
                symbolcls = self.algebra.codegen_symbolcls
                if hasattr(symbolcls, "__self__"):
                    symbolcls = symbolcls.__self__
                if isinstance(v, symbolcls):
                    return True
        
        # No symbolic values found
        return False
    
    
    @property
    def free_symbols(self) -> Set:
        """
        Get the set of free symbols present in the multivector.
        
        Returns:
            Set of sympy Symbols
        """
        symbols = set()
        for v in self.values():
            if hasattr(v, "free_symbols"):
                symbols |= v.free_symbols
        return symbols
    
    @property
    def type_number(self) -> int:
        """
        Compute a unique type number based on present basis blades.
        
        Returns:
            Integer representing the type of multivector
        """
        bitstr = "".join(
            "1" if i in self.keys() else "0" 
            for i in reversed(self.algebra.canon2bin.values())
        )
        return int(bitstr, 2)
    
    # Custom attribute access for basis blades
    def __getattr__(self, name: str) -> Any:
        """
        Allow direct access to basis blade coefficients as attributes.
        
        Args:
            name: Name of attribute to access
            
        Returns:
            Value of the basis blade coefficient
            
        Raises:
            AttributeError: If the attribute is not a valid basis blade
        """
        # Handle special case for scalar
        if name == 'e':
            scalar_idx = self._keys.index(0) if 0 in self._keys else -1
            return self._values[scalar_idx] if scalar_idx >= 0 else 0
            
        # Try to get from basis blades
        if hasattr(self.algebra, 'canon2bin'):
            blade_key = self.algebra.canon2bin.get(name)
            if blade_key is not None:
                return self._values[self._keys.index(blade_key)] if blade_key in self._keys else 0
        
        raise AttributeError(f"'{type(self).__name__}' object has no attribute '{name}'")
    
    # === Core Methods ===
    
    def keys(self) -> Tuple[int, ...]:
        """
        Get the binary keys for basis blades in this multivector.
        
        Returns:
            Tuple of integer keys representing basis blades
        """
        return self._keys
    
    def values(self) -> Union[List[Any], np.ndarray]:
        """
        Get the coefficients for basis blades in this multivector.
        
        Returns:
            List or array of coefficients
        """
        return self._values
    
    
    def itermv(self) -> Generator['MultiVector', None, None]:
        """
        Iterate over the multivector elements if it has a shape.
        
        Yields:
            MultiVector elements for each index in the shape
        """
        if not hasattr(self, 'shape') or not self.shape:
            yield self
            return
        
        # Get the shape of the multivector
        shape = self.shape
        
        # Iterate over all indices in the shape
        if isinstance(self._values, np.ndarray) and self._values.ndim > 1:
            for idx in np.ndindex(shape):
                # Extract the values for this index
                values = self._values[(slice(None),) + idx]
                yield type(self).fromkeysvalues(self.algebra, self._keys, values)
        else:
            # Handle non-numpy values
            for i in range(shape[0] if shape else 1):
                # Extract the values for this index
                values = [v[i] if isinstance(v, (list, tuple)) else v for v in self._values]
                yield type(self).fromkeysvalues(self.algebra, self._keys, values)
                
    def map(self, func: Callable) -> 'MultiVector':
        """
        Apply a function to each component of a multivector.
        
        Args:
            func: Function to apply to each component value or (key, value) pair
            
        Returns:
            A new multivector with transformed values
        """
        if callable(func) and func.__code__.co_argcount == 2:
            # Function takes key and value
            new_values = [func(k, v) for k, v in zip(self._keys, self._values)]
        else:
            # Function takes only value
            new_values = [func(v) for v in self._values]
        
        # Create a new multivector with the transformed values
        return type(self).fromkeysvalues(self.algebra, self._keys, new_values)

    def filter(self, pred: Callable) -> 'MultiVector':
        """
        Filter components of a multivector based on a predicate function.
        
        Args:
            pred: Predicate function that takes a value or (key, value) and returns True/False
            
        Returns:
            A new multivector with only the components for which pred returns True
        """
        if callable(pred) and pred.__code__.co_argcount == 2:
            # Predicate takes key and value
            filtered_items = [(k, v) for k, v in zip(self._keys, self._values) if pred(k, v)]
        else:
            # Predicate takes only value
            filtered_items = [(k, v) for k, v in zip(self._keys, self._values) if pred(v)]
        
        if not filtered_items:
            return type(self).fromkeysvalues(self.algebra, (), [])
            
        keys, values = zip(*filtered_items) if filtered_items else ((), [])
        return type(self).fromkeysvalues(self.algebra, keys, list(values))        
    
    def exp(self) -> 'MultiVector':
        """
        Compute the exponential of a multivector.
        
        Returns:
            A new multivector representing exp(self)
        """
        # Special case for zero
        if not self._keys or all(v == 0 for v in self._values):
            return self.algebra.scalar(1)
            
        # For vectors and bivectors, use specialized formulas
        if len(self.grades) == 1 and self.grades[0] in (1, 2):
            # Square the multivector
            mv2 = self * self
            
            # Check if it's nilpotent (x² = 0)
            if mv2.grades == (0,) and abs(mv2.e) < 1e-10:
                return self.algebra.scalar(1) + self
                
            # Check if it has a negative square (typical for bivectors)
            if mv2.grades == (0,) and mv2.e < 0:
                # Use formula exp(x) = cos(|x|) + x/|x| * sin(|x|)
                import math
                alpha = math.sqrt(-mv2.e)
                return self.algebra.scalar(math.cos(alpha)) + (self / alpha) * math.sin(alpha)
                
            # If it has a positive square, use hyperbolic functions
            if mv2.grades == (0,) and mv2.e > 0:
                import math
                alpha = math.sqrt(mv2.e)
                return self.algebra.scalar(math.cosh(alpha)) + (self / alpha) * math.sinh(alpha)
                
        # For general multivectors, use Taylor series
        result = self.algebra.scalar(1)
        term = self.algebra.scalar(1)
        factorial = 1
        
        # Compute the series up to 15 terms
        for i in range(1, 16):
            factorial *= i
            term = term * self
            result += term / factorial
            
            # Stop if the term becomes negligible
            if hasattr(term, '_values') and all(abs(v) < 1e-10 for v in term._values):
                break
                
        return result          
                
    
    def items(self) -> Generator[Tuple[int, Any], None, None]:
        """
        Get iterator over (key, value) pairs.
        
        Yields:
            Pairs of (key, coefficient)
        """
        yield from zip(self._keys, self._values)
    
    def copy(self) -> MultiVector:
        """
        Create a copy of this multivector.
        
        Returns:
            A new multivector with the same content
        """
        values_copy = self._values.copy() if isinstance(self._values, np.ndarray) else list(self._values)
        return self.__class__.fromkeysvalues(self.algebra, self._keys, values_copy)
    
    def __len__(self) -> int:
        """
        Get the number of components in the multivector.
        
        Returns:
            Number of basis blades with coefficients
        """
        return len(self._values)
    
    # === Grade Filtering ===
    
    def grade(self, *grades: int) -> 'MultiVector':
        """
        Extract specific grades from the multivector.
        
        Args:
            *grades: Grades to extract (0=scalar, 1=vector, 2=bivector, etc.)
            
        Returns:
            A new multivector containing only the specified grades
        """
        if not grades:
            return type(self).fromkeysvalues(self.algebra, (), [])
        
        # Handle case where grades is a tuple within a tuple
        if len(grades) == 1 and isinstance(grades[0], tuple):
            grades = grades[0]
        
        # Find keys of specified grades
        matching_keys = [k for k in self._keys if bin(k).count('1') in grades]
        if not matching_keys:
            return type(self).fromkeysvalues(self.algebra, (), [])
            
        # Extract values for matching keys
        matching_indices = [self._keys.index(k) for k in matching_keys]
        
        # Handle both numpy arrays and regular lists
        if isinstance(self._values, np.ndarray):
            matching_values = self._values[matching_indices]
        else:
            matching_values = [self._values[i] for i in matching_indices]
        
        return type(self).fromkeysvalues(self.algebra, tuple(matching_keys), matching_values)
    
    def scalar(self) -> MultiVector:
        """
        Extract the scalar part (grade 0).
        
        Returns:
            A new multivector containing only the scalar part
        """
        return self.grade(0)
    
    def vector(self) -> MultiVector:
        """
        Extract the vector part (grade 1).
        
        Returns:
            A new multivector containing only the vector part
        """
        return self.grade(1)
    
    def bivector(self) -> MultiVector:
        """
        Extract the bivector part (grade 2).
        
        Returns:
            A new multivector containing only the bivector part
        """
        return self.grade(2)
    
    # === Geometric Operations ===
    
    def reverse(self) -> MultiVector:
        """
        Compute the reverse of this multivector.
        
        The reverse flips the order of basis vectors, changing signs for grades 2 and 3 (mod 4).
        
        Returns:
            A new multivector representing the reverse
        """
        result = {}
        for k, v in self.items():
            grade = bin(k).count('1')  # Calculate the grade
            sign = -1 if grade % 4 in (2, 3) else 1  # Grades 2 and 3 (mod 4) change sign
            result[k] = -v if sign < 0 else v  # Apply sign change
            
        # Convert result to a proper multivector
        keys, values = zip(*result.items()) if result else ((), [])
        return type(self).fromkeysvalues(self.algebra, keys, list(values))
    
    def involute(self) -> MultiVector:
        """
        Compute the grade involution of this multivector.
        
        The grade involution negates odd-grade components.
        
        Returns:
            A new multivector representing the grade involution
        """
        return self.algebra.involute(self)
    
    def conjugate(self) -> MultiVector:
        """
        Compute the Clifford conjugate of this multivector.
        
        The Clifford conjugate combines reversion and grade involution.
        
        Returns:
            A new multivector representing the Clifford conjugate
        """
        return self.algebra.conjugate(self)
    
    def dual(self) -> MultiVector:
        """
        Compute the dual of this multivector.
        
        The dual maps a multivector to its orthogonal complement using the pseudoscalar.
        
        Returns:
            A new multivector representing the dual
        """
        return self.algebra.dual(self)
    
    def inverse(self) -> 'MultiVector':
        """
        Compute the inverse of this multivector.
        
        Returns:
            A new multivector representing the inverse
            
        Raises:
            ZeroDivisionError: If the multivector is not invertible
        """
        return self.algebra.inv(self)

    def inv(self) -> 'MultiVector':
        """
        Compute the inverse of this multivector.

        The inverse of a multivector a is the multivector a^-1 such that a * a^-1 = a^-1 * a = 1.

        Returns:
            A new multivector representing the inverse
        
        Raises:
            ZeroDivisionError: If the multivector has no valid inverse
        """
        # Direct scalar inverse for efficiency and robustness
        if hasattr(self, '_keys') and isinstance(self._keys, (tuple, list)) and len(self._keys) == 1 and self._keys[0] == 0:
            value = self._values[0]
            if value == 0:
                raise ZeroDivisionError("Cannot compute inverse of zero scalar")
            return type(self).fromkeysvalues(self.algebra, (0,), [1.0 / value])
        
        # For everything else, use the algebra's inverse operation
        try:
            return self.algebra.inv(self)
        except (TypeError, ValueError) as e:
            # Special case for single-grade multivectors
            if self.grades == (0,):  # Pure scalar
                if len(self._values) == 1 and self._values[0] != 0:
                    return type(self).fromkeysvalues(self.algebra, (0,), [1.0 / self._values[0]])
            elif self.grades == (1,):  # Pure vector
                # For vectors, compute |v|^2
                norm_sq = sum(v**2 for v in self._values)
                if norm_sq == 0:
                    raise ZeroDivisionError("Vector has zero norm")
                return self / norm_sq
            
            # If we couldn't handle it with special cases, re-raise
            raise
    
        
    def norm(self) -> float:
        """
        Compute the norm (magnitude) of this multivector.
        
        The norm is computed as sqrt(|a * ~a|) where ~a is the reverse.
        For vectors this is the usual Euclidean norm.
        
        Returns:
            The norm as a float value
        """
        # Calculate |A * ~A| = |A|²
        norm_sq = (self * self.reverse()).scalar()
        if not norm_sq or len(norm_sq._values) == 0:
            return 0.0
            
        # Extract the scalar value
        if hasattr(norm_sq, 'e'):
            norm_sq_val = norm_sq.e
        else:
            # If e attribute is not available, try to access the first value
            norm_sq_val = norm_sq._values[0] if len(norm_sq._values) > 0 else 0
            
        # Take the square root
        import math
        if norm_sq_val < 0:
            # Handle negative squared norm (can happen with certain signatures)
            return math.sqrt(abs(norm_sq_val))
        else:
            return math.sqrt(norm_sq_val)

    def normalized(self) -> MultiVector:
        """
        Return a normalized version of this multivector.
        
        The normalization divides the multivector by its norm.
        
        Returns:
            A new multivector with unit norm
            
        Raises:
            ZeroDivisionError: If the multivector has zero norm
        """
        norm_value = self.norm()
        if norm_value == 0:
            raise ZeroDivisionError("Cannot normalize a multivector with zero norm")
        return self / norm_value
    
        
    # === Arithmetic Operations ===
    
    def __eq__(self, other: Any) -> bool:
        """
        Check if this multivector equals another multivector or scalar.
        
        Args:
            other: Another multivector or scalar value
            
        Returns:
            True if equal, False otherwise
        """
        if isinstance(other, (int, float, complex)):
            # For scalar comparison, check if this is a scalar with the right value
            if not self._keys:
                return other == 0
            return len(self._keys) == 1 and self._keys[0] == 0 and self._values[0] == other
            
        if isinstance(other, MultiVector):
            # For multivector comparison, compare all components
            if len(self._keys) != len(other._keys):
                return False
                
            self_dict = dict(zip(self._keys, self._values))
            other_dict = dict(zip(other._keys, other._values))
            
            # Check each key in either multivector
            all_keys = set(self_dict) | set(other_dict)
            return all(self_dict.get(k, 0) == other_dict.get(k, 0) for k in all_keys)
            
        return NotImplemented
    
    def __add__(self, other: Union[MultiVector, int, float, complex]) -> MultiVector:
        """
        Add this multivector to another multivector or scalar.
        
        Args:
            other: Another multivector or scalar value
            
        Returns:
            A new multivector representing the sum
        """
        if isinstance(other, (int, float, complex)):
            # For scalar addition, add to the scalar part or create it
            result_dict = dict(zip(self._keys, self._values))
            result_dict[0] = result_dict.get(0, 0) + other
            
            # Filter out zero values
            result_dict = {k: v for k, v in result_dict.items() if v != 0}
            keys, values = zip(*result_dict.items()) if result_dict else ((), [])
            
            return type(self).fromkeysvalues(self.algebra, keys, list(values))
            
        if isinstance(other, MultiVector):
            # For multivector addition, combine all components
            if self.algebra != other.algebra:
                raise ValueError("Cannot add multivectors from different algebras")
                
            result_dict = dict(zip(self._keys, self._values))
            for k, v in zip(other._keys, other._values):
                result_dict[k] = result_dict.get(k, 0) + v
                
            # Filter out zero values
            result_dict = {k: v for k, v in result_dict.items() if v != 0}
            keys, values = zip(*result_dict.items()) if result_dict else ((), [])
            
            return type(self).fromkeysvalues(self.algebra, keys, list(values))
            
        return NotImplemented
    
    def __radd__(self, other: Union[int, float, complex]) -> MultiVector:
        """
        Add a scalar to this multivector (right addition).
        
        Args:
            other: A scalar value
            
        Returns:
            A new multivector representing the sum
        """
        return self + other
    
    def __sub__(self, other: Union[MultiVector, int, float, complex]) -> MultiVector:
        """
        Subtract another multivector or scalar from this multivector.
        
        Args:
            other: Another multivector or scalar value
            
        Returns:
            A new multivector representing the difference
        """
        if isinstance(other, (int, float, complex)):
            # For scalar subtraction, subtract from the scalar part or create it
            result_dict = dict(zip(self._keys, self._values))
            result_dict[0] = result_dict.get(0, 0) - other
            
            # Filter out zero values
            result_dict = {k: v for k, v in result_dict.items() if v != 0}
            keys, values = zip(*result_dict.items()) if result_dict else ((), [])
            
            return type(self).fromkeysvalues(self.algebra, keys, list(values))
            
        if isinstance(other, MultiVector):
            # For multivector subtraction, combine all components
            if self.algebra != other.algebra:
                raise ValueError("Cannot subtract multivectors from different algebras")
                
            result_dict = dict(zip(self._keys, self._values))
            for k, v in zip(other._keys, other._values):
                result_dict[k] = result_dict.get(k, 0) - v
                
            # Filter out zero values
            result_dict = {k: v for k, v in result_dict.items() if v != 0}
            keys, values = zip(*result_dict.items()) if result_dict else ((), [])
            
            return type(self).fromkeysvalues(self.algebra, keys, list(values))
            
        return NotImplemented
    
    def __rsub__(self, other: Union[int, float, complex]) -> MultiVector:
        """
        Subtract this multivector from a scalar (right subtraction).
        
        Args:
            other: A scalar value
            
        Returns:
            A new multivector representing the difference
        """
        # Create a negated copy, then add the scalar
        neg_values = [-v for v in self._values]
        result_dict = dict(zip(self._keys, neg_values))
        result_dict[0] = result_dict.get(0, 0) + other
        
        # Filter out zero values
        result_dict = {k: v for k, v in result_dict.items() if v != 0}
        keys, values = zip(*result_dict.items()) if result_dict else ((), [])
        
        return type(self).fromkeysvalues(self.algebra, keys, list(values))
    
    def __mul__(self, other: Union[MultiVector, int, float, complex]) -> MultiVector:
        """
        Multiply this multivector by another multivector or scalar.
        
        Args:
            other: Another multivector or scalar value
            
        Returns:
            A new multivector representing the geometric product
        """
        if isinstance(other, (int, float, complex)):
            # For scalar multiplication, multiply all components
            new_values = [v * other for v in self._values]
            return type(self).fromkeysvalues(self.algebra, self._keys, new_values)
            
        if isinstance(other, MultiVector):
            # For multivector multiplication, use the algebra's geometric product
            result = self.algebra.gp(self, other)
            
            # Special case for vector * vector = scalar
            # Check if both are pure vectors (grade 1) and identical
             # Special case for pure-vector geometric product
            if self.grades == (1,) and other.grades == (1,):
                k1, k2 = self._keys[0], other._keys[0]
                key = k1 ^ k2
                # Use sign table (0 if orthogonal → zero bivector)
                sign = self.algebra.signs.get((k1, k2), 0)
                return type(self).fromkeysvalues(self.algebra, (key,), [sign])

            
            # Special case: if we got an empty result, return a zero scalar
            if not result._keys or len(result._keys) == 0:
                return type(self).fromkeysvalues(self.algebra, (0,), [0])
                
            return result
            
        return NotImplemented
    
    def __rmul__(self, other: Union[int, float, complex]) -> MultiVector:
        """
        Multiply a scalar by this multivector (right multiplication).
        
        Args:
            other: A scalar value
            
        Returns:
            A new multivector representing the product
        """
        return self * other
    
    def __truediv__(self, other: Union[MultiVector, int, float, complex]) -> MultiVector:
        """
        Divide this multivector by another multivector or scalar.
        
        Args:
            other: Another multivector or scalar value
            
        Returns:
            A new multivector representing the quotient
            
        Raises:
            ZeroDivisionError: If dividing by zero
        """
        if isinstance(other, (int, float, complex)):
            # For scalar division, divide all components
            if other == 0:
                raise ZeroDivisionError("Division by zero")
                
            new_values = [v / other for v in self._values]
            return type(self).fromkeysvalues(self.algebra, self._keys, new_values)
            
        if isinstance(other, MultiVector):
            # For multivector division, use the algebra's division operation
            return self.algebra.div(self, other)
            
        return NotImplemented
    
    def __neg__(self) -> MultiVector:
        """
        Negate this multivector.
        
        Returns:
            A new multivector with all components negated
        """
        new_values = [-v for v in self._values]
        return type(self).fromkeysvalues(self.algebra, self._keys, new_values)
    
    def __invert__(self) -> MultiVector:
        """
        Compute the reverse of this multivector (~ operator).
        
        Returns:
            A new multivector representing the reverse
        """
        return self.reverse()
    
    def __pow__(self, power: int) -> 'MultiVector':
        """
        Raise this multivector to a power.
        
        Args:
            power: Integer exponent
            
        Returns:
            A new multivector representing the power
            
        Raises:
            ValueError: If the power is negative
        """
        if isinstance(power, int):
            if power < 0:
                # For negative powers, compute the inverse and then positive power
                return self.inv() ** (-power)
            
            if power == 0:
                return self.algebra.scalar(1)
                
            if power == 1:
                return self.copy()
                
            # Use binary exponentiation for efficiency
            result = self.algebra.scalar(1)
            base = self.copy()
            
            while power > 0:
                if power & 1:
                    result = result * base
                base = base * base
                power >>= 1
                
            return result
        elif power == 0.5:
            # Square root
            return self.sqrt()
        elif power == -0.5:
            # Inverse square root
            return self.inv().sqrt()
        else:
            raise ValueError(f"Unsupported power: {power}")
        
        
    def asfullmv(self, canonical: bool = True) -> 'MultiVector':
        """
        Convert to a full multivector with all possible components.
        
        Args:
            canonical: If True, use canonical ordering of components
            
        Returns:
            A new multivector with all possible components
        """
        # Create a dictionary with all components (initially zero)
        if canonical:
            # Use canonical ordering from algebra
            values = np.zeros(len(self.algebra))
            
            # Fill in the non-zero components
            for idx, k in enumerate(self.algebra.canon2bin.values()):
                if k in self._keys:
                    values[idx] = self._values[self._keys.index(k)]
                    
            return type(self).fromkeysvalues(self.algebra, tuple(self.algebra.canon2bin.values()), values)
        else:
            # Use binary ordering
            values = np.zeros(len(self.algebra))
            
            # Fill in the non-zero components
            for k, v in zip(self._keys, self._values):
                if k in set(range(2**self.algebra.d)):
                    values[k] = v
                    
            return type(self).fromkeysvalues(self.algebra, tuple(range(2**self.algebra.d)), values)       
    
    # === Geometric Algebra Specific Operations ===
    
            
    def __xor__(self, other: MultiVector) -> MultiVector:
        if isinstance(other, MultiVector) and self.grades == (1,) and other.grades == (1,):
            key = self._keys[0] ^ other._keys[0]
            return type(self).fromkeysvalues(self.algebra, (key,), [1])
        return self.algebra.op(self, other)
    
    
 
    
    
    def __or__(self, other: MultiVector) -> MultiVector:
        """
        Compute the inner product with another multivector.
        
        Args:
            other: Another multivector
            
        Returns:
            A new multivector representing the inner product
        """
        if isinstance(other, MultiVector):
            return self.algebra.ip(self, other)
        return NotImplemented
    
    def __and__(self, other: MultiVector) -> MultiVector:
        """
        Compute the regressive product with another multivector.
        
        Args:
            other: Another multivector
            
        Returns:
            A new multivector representing the regressive product
        """
        if isinstance(other, MultiVector):
            return self.algebra.rp(self, other)
        return NotImplemented
    
    def __rshift__(self, other: MultiVector) -> MultiVector:
        """
        Compute the sandwich product with another multivector.
        
        The sandwich product is x >> y = x * y * ~x, frequently used for rotations.
        
        Args:
            other: Another multivector
            
        Returns:
            A new multivector representing the sandwich product
        """
        if isinstance(other, MultiVector):
            return self.algebra.sw(self, other)
        return NotImplemented
    
    def __matmul__(self, other: MultiVector) -> MultiVector:
        """
        Compute the projection of another multivector onto this multivector.
        
        Args:
            other: Another multivector
            
        Returns:
            A new multivector representing the projection
        """
        if isinstance(other, MultiVector):
            return self.algebra.proj(self, other)
        return NotImplemented
    
    # === Indexing and Array Interface ===
    
    def __getitem__(self, item: Union[int, Tuple, slice]) -> 'MultiVector':
        """
        Allow flexible indexing into the multivector's coefficients.
        
        Supports:
        - Single integer index
        - Tuple indices
        - Slicing
        - Extracting array-valued multivectors
        
        Args:
            item: Index or slice
            
        Returns:
            A new multivector or component value
        """
        # Convert single index to tuple
        if not isinstance(item, tuple):
            item = (item,)
        
        vals = self.values()
        
        # Handle array-like values
        if isinstance(vals, (list, tuple)):
            # For single integer index, return that component
            if len(item) == 1 and isinstance(item[0], int):
                if 0 <= item[0] < len(self._keys):
                    return type(self).fromkeysvalues(
                        self.algebra, 
                        (self._keys[item[0]],), 
                        [vals[item[0]]]
                    )
                else:
                    raise IndexError(f"Index {item[0]} out of range for multivector with {len(self._keys)} components")
            
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
                try:
                    new_vals = vals[item]
                    # Important: Preserve shape information for sliced arrays
                    if isinstance(new_vals, np.ndarray):
                        # Create result multivector
                        mv = type(self).fromkeysvalues(
                            self.algebra, 
                            self._keys, 
                            new_vals
                        )
                        # Explicitly set the shape attribute
                        mv.shape = (new_vals.size,) if new_vals.ndim == 1 else new_vals.shape
                        return mv
                except (IndexError, TypeError):
                    # Fallback for invalid indexing
                    new_vals = vals
            else:
                # For 1D arrays, handle integer indexing
                if len(item) == 1 and isinstance(item[0], int):
                    if 0 <= item[0] < len(self._keys):
                        return type(self).fromkeysvalues(
                            self.algebra, 
                            (self._keys[item[0]],), 
                            [vals[item[0]]]
                        )
                    else:
                        raise IndexError(f"Index {item[0]} out of range for multivector with {len(self._keys)} components")
                try:
                    new_vals = vals[item]
                    if isinstance(new_vals, np.ndarray) and new_vals.ndim == 1:
                        # Ensure proper shape for 1D array result
                        mv = type(self).fromkeysvalues(self.algebra, self._keys, new_vals)
                        mv.shape = (new_vals.size,)  # Explicitly set shape for 1D arrays
                        return mv
                except (IndexError, TypeError):
                    new_vals = vals
        else:
            # For scalar values, just return them
            new_vals = vals
        
        # Create result multivector
        return type(self).fromkeysvalues(self.algebra, self._keys, new_vals)
    
    def __setitem__(self, key: Union[int, Tuple, slice], value: Any) -> None:
        """
        Handle assignment to MultiVector objects with improved NumPy interoperability.
        
        Supports:
        - Assigning values to specific components
        - Assigning values to array-valued multivectors
        - Assigning values from NumPy arrays
        
        Args:
            key: Index or slice
            value: Value to assign
        """
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
                idx = key[0] if len(key) == 1 else key
                if isinstance(idx, int) and 0 <= idx < len(self._values):
                    self._values[idx] = value
                elif isinstance(idx, slice):
                    indices = range(*idx.indices(len(self._values)))
                    values = [value] * len(indices) if not isinstance(value, (list, tuple, np.ndarray)) else value
                    for i, v in zip(indices, values):
                        self._values[i] = v
            else:
                # Fallback for other unexpected cases
                self._values[key[0] if len(key) == 1 else key] = value
    
    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        """
        Support for NumPy's universal functions.
        
        This method enables multivectors to work with NumPy's mathematical functions.
        
        Args:
            ufunc: NumPy universal function
            method: Method being called ('__call__', etc.)
            *inputs: Input arguments
            **kwargs: Keyword arguments
            
        Returns:
            Result of applying the ufunc, wrapped in a multivector if appropriate
        """
        # Extract the multivector values and apply the ufunc
        if method == '__call__':
            # Extract values from multivector inputs
            args = []
            for input_arg in inputs:
                if isinstance(input_arg, MultiVector):
                    args.append(input_arg.values())
                else:
                    args.append(input_arg)
                    
            # Apply the ufunc to the values
            result = ufunc(*args, **kwargs)
            
            # Wrap the result in a new multivector
            return type(self).fromkeysvalues(self.algebra, self._keys, result)
            
        return NotImplemented
    
        
    # === Representation and String Methods ===
    
    def __str__(self) -> str:
        """
        Convert to a human-readable string representation.
        
        Returns:
            String representation of the multivector
        """
        if not self._values or len(self._values) == 0:
            return "0"
            
        def format_val(val: Any) -> str:
            """Format a value for display"""
            s = str(val)
            if isinstance(val, sympy.Expr):
                return s if val.is_Symbol else f"({s})"
            if isinstance(val, float):
                return f"{val:.3g}"
            return s
        
        # Make sure we have the proper mapping attribute
        if not hasattr(self.algebra, '_bin2canon_prettystr'):
            # If missing, create a simple fallback
            self.algebra._bin2canon_prettystr = {}
            for k in self._keys:
                if k == 0:
                    self.algebra._bin2canon_prettystr[k] = "1"
                else:
                    indices = [i for i in range(self.algebra.d) if (k & (1 << i))]
                    idx_str = ''.join(str(i+1) for i in indices)
                    self.algebra._bin2canon_prettystr[k] = f"e{idx_str}"
        
        # Create pretty string representation
        pretty_vals = {self.algebra._bin2canon_prettystr.get(k, str(k)): v for k, v in self.items()}
        components = []
        
        for blade, val in pretty_vals.items():
            # Check if the value is non-zero
            if hasattr(val, "any"):  # For NumPy arrays
                keep = val.any()
            else:
                keep = bool(val)
            if not keep:
                continue
                
            if blade == "1":
                components.append(format_val(val))
            else:
                if val == 1:
                    components.append(blade)
                elif val == -1:
                    components.append(f"-{blade}")
                else:
                    components.append(f"{format_val(val)}{blade}")
        
        return " + ".join(components).replace("+ -", "- ") if components else "0"
    
    def __repr__(self) -> str:
        """
        Generate a detailed string representation for debugging.
        
        Returns:
            Detailed string representation
        """
        return f"MultiVector({self.algebra}, {dict(self.items())})"
    
    # === Transformation Methods ===
    
    def dual(self, kind: str = 'auto') -> MultiVector:
        """
        Compute the dual of this multivector.
        
        The dual maps a multivector to its orthogonal complement using the pseudoscalar.
        
        Args:
            kind: Type of dual ('polarity', 'hodge', or 'auto')
            
        Returns:
            A new multivector representing the dual
        """
        if kind == 'polarity' or (kind == 'auto' and self.algebra.r == 0):
            return self.algebra.polarity(self)
        elif kind == 'hodge' or (kind == 'auto' and self.algebra.r == 1):
            return self.algebra.hodge(self)
        elif kind == 'auto':
            raise ValueError('Cannot select a suitable dual in auto mode for this algebra.')
        else:
            raise ValueError(f'No dual found for kind={kind}.')
    
    def undual(self, kind: str = 'auto') -> MultiVector:
        """
        Compute the undual (inverse of dual) of this multivector.
        
        Args:
            kind: Type of undual ('polarity', 'hodge', or 'auto')
            
        Returns:
            A new multivector representing the undual
        """
        if kind == 'polarity' or (kind == 'auto' and self.algebra.r == 0):
            return self.algebra.unpolarity(self)
        elif kind == 'hodge' or (kind == 'auto' and self.algebra.r == 1):
            return self.algebra.unhodge(self)
        elif kind == 'auto':
            raise ValueError('Cannot select a suitable undual in auto mode for this algebra.')
        else:
            raise ValueError(f'No undual found for kind={kind}.')
    
    def tosympy(self) -> Any:
        """
        Convert to a sympy expression if possible.
        
        Returns:
            A sympy expression representing this multivector
        """
        if len(self._keys) == 1 and self._keys[0] == 0:
            # For a scalar, just return the value
            val = self._values[0]
            return val if isinstance(val, sympy.Expr) else sympy.S(val)
        else:
            # For non-scalar, return a dict representation
            return {k: (v if isinstance(v, sympy.Expr) else sympy.S(v)) for k, v in self.items()}
        
    def __deepcopy__(self, memo):
        from copy import deepcopy
        keys = tuple(self._keys)
        values = deepcopy(list(self._values), memo)
        return type(self).fromkeysvalues(self.algebra, keys, values)
        
    def __deepcopy__(self, memo):
        from copy import deepcopy
        keys = tuple(self._keys)
        values = deepcopy(list(self._values), memo)
        return type(self).fromkeysvalues(self.algebra, keys, values)


# Define aliases for backward compatibility
MultiVectorBase = MultiVector
MultiVectorCreation = MultiVector
MultiVectorIndexing = MultiVector
MultiVectorOperations = MultiVector
MultiVectorTransforms = MultiVector

# Export the MultiVector class
__all__ = ['MultiVector']
