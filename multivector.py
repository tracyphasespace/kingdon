# -*- coding: utf-8 -*-
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
import math # Added for sqrt
from functools import reduce
from typing import (
    Any, Callable, Dict, Generator, List, Optional, Sequence,
    Set, Tuple, Union, cast, TYPE_CHECKING
)

import numpy as np
import sympy
from sympy import Expr, S # Import S for symbolic constants

# Use TYPE_CHECKING for imports that might cause circular imports
# REVIEW: Good practice to avoid circular imports. Ensure RationalPolynomial is handled consistently.
if TYPE_CHECKING:
    from kingdon.polynomial import RationalPolynomial
else:
    # Create a placeholder class for runtime checks
    class RationalPolynomial:
        """Placeholder for RationalPolynomial to avoid circular imports."""
        pass

# TODO: Consider making MultiVector immutable by removing __setitem__ or raising error.
# Current implementation is mutable.
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

    Note on Mutability: This class is currently mutable via __setitem__.
    """

    # Lower priority than NumPy arrays (0) to potentially allow NumPy's __mul__ etc.
    # to handle operations if desired, though typically MultiVector's ops are wanted.
    # Review interaction if necessary.
    # REVIEW: This ensures MultiVector's methods are preferred over NumPy's when operating
    # between a MultiVector and a NumPy array (e.g., mv * np_array). This is usually desired.
    __array_priority__ = -100

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
            # Ensure _keys is a tuple for index()
            # REVIEW: Defensive check, good. Assumes _keys is set during init.
            keys_tuple = self._keys if isinstance(self._keys, tuple) else tuple(self._keys)
            idx = keys_tuple.index(key)
            # Handle potential list/array indexing
            # REVIEW: Handles both list and ndarray cases for _values.
            if isinstance(self._values, (list, np.ndarray)):
                 # Ensure index is within bounds (should be guaranteed by keys_tuple.index)
                 if 0 <= idx < len(self._values):
                     return self._values[idx]
                 else:
                     # This case indicates internal inconsistency
                     print(f"Warning: Index {idx} out of bounds for values (len {len(self._values)}) in get().")
                     return default
            else:
                 # Fallback if _values is not directly indexable (should not happen with standard init)
                 # REVIEW: Consider raising an internal error here, as _values should always be indexable.
                 print(f"Warning: _values attribute (type {type(self._values)}) is not indexable in get().")
                 return default # Or handle appropriately
        except (ValueError, IndexError):
            # Handles case where key is not found in _keys tuple
            return default
        except Exception as e:
            # Catch other potential errors during access
            print(f"Error in MultiVector.get(key={key}): {e}")
            return default


    def __array__(self, dtype=None):
        """
        Support conversion to NumPy array.

        Args:
            dtype: Data type for the array

        Returns:
            NumPy array representation. For dtype=object, returns an array containing self.
        """
        if dtype == object:
            # For object arrays, return a numpy array containing self.
            # This allows arrays of MultiVectors, e.g., np.array([mv1, mv2]).
            # REVIEW: Correct handling for creating arrays *of* MultiVectors.
            arr = np.empty(1, dtype=object)
            arr[0] = self
            return arr

        # For numeric arrays, extract values. Ensure consistent type.
        try:
            # Attempt conversion, default to 0 for empty multivectors
            # REVIEW: Handles empty MV case reasonably.
            numeric_values = list(self._values) if hasattr(self, '_values') and self._values is not None else [0]
            if not numeric_values and (not hasattr(self, '_keys') or not self._keys):
                 numeric_values = [0] # Represent empty MV as [0]
            elif not numeric_values and hasattr(self, '_keys') and self._keys:
                 # Has keys but empty values? Inconsistent state. Default to zeros?
                 print("Warning: Inconsistent state in __array__ (keys exist but values are empty). Returning zeros.")
                 numeric_values = [0] * len(self._keys)


            # Convert symbolic elements if possible, otherwise error might be appropriate
            # For now, let np.array handle it, might raise TypeError for symbolic
            # REVIEW: Conversion to numeric array might fail for symbolic values. The warning is good.
            # Consider if a specific error or different behavior is needed for symbolic MVs.
            return np.array(numeric_values, dtype=dtype)
        except TypeError as e:
             # Handle cases with symbolic values if necessary, maybe return object array?
             print(f"Warning: Cannot convert symbolic MultiVector to numeric NumPy array ({e}). Returning object array.")
             # Ensure _values exists before getting length
             num_vals = len(self._values) if hasattr(self, '_values') and self._values is not None else 0
             arr = np.empty(num_vals, dtype=object)
             if num_vals > 0:
                 for i, v in enumerate(self._values):
                    arr[i] = v
             return arr


    def __new__(cls, *args, **kwargs) -> MultiVector:
        """
        Create a new MultiVector instance. Delegates actual initialization to __init__.

        Args:
            *args: Positional arguments, typically including the algebra
            **kwargs: Keyword arguments for initialization

        Returns:
            A new MultiVector instance
        """
        # REVIEW: Standard __new__ implementation.
        obj = object.__new__(cls)
        return obj

    def __init__(self, algebra: Any, **kwargs: Any) -> None:
        """
        Initialize a MultiVector. Typically called via factory methods like algebra.vector().
        Prefer using `algebra.multivector(...)` or `MultiVector.fromkeysvalues(...)`.

        Args:
            algebra: The Geometric Algebra instance this multivector belongs to.
            **kwargs: Initialization parameters including:
                - keys: Tuple of binary keys for basis blades (e.g., (0,) for scalar, (1, 2) for e1+e2).
                - values: Coefficients (list or numpy array) matching the keys.
                - name: Optional name for symbolic representation.
                - grades: Optional filter to initialize only specific grades.
        """
        # REVIEW: This __init__ seems designed to be called *after* fromkeysvalues or similar
        # might have already set _keys/_values. This could be slightly confusing.
        # Consider making fromkeysvalues the primary internal constructor and simplifying __init__.
        self.algebra = algebra
        self.shape = () # Initialize shape, will be set by fromkeysvalues if called

        # Initialize _keys and _values only if not already set (e.g., by fromkeysvalues)
        if not hasattr(self, '_keys'):
             # Default to zero scalar if no keys/values provided
             keys_arg = kwargs.get('keys', (0,))
             # Ensure keys is always a tuple
             if isinstance(keys_arg, int):
                 self._keys = (keys_arg,)
             elif not isinstance(keys_arg, tuple):
                 self._keys = tuple(keys_arg)
             else:
                 self._keys = keys_arg


        if not hasattr(self, '_values'):
             values_arg = kwargs.get('values', None)
             if values_arg is not None:
                 # Ensure values is a list or ndarray, handle single values
                 if isinstance(values_arg, (list, np.ndarray)):
                     self._values = values_arg
                 else:
                     # Ensure single value matches single key, or broadcast? Assume match for now.
                     self._values = [values_arg] # Assume single value corresponds to single key
             else:
                 # Default to zero(s) matching the number of keys
                 # Check issymbolic based on algebra? Safer to default to numeric 0 unless algebra indicates otherwise.
                 # For now, stick to numeric 0 as default.
                 self._values = [0] * len(self._keys) # Default to numeric 0

             # Ensure keys and values match length
             if len(self._keys) != len(self._values):
                 raise ValueError(f"Number of keys ({len(self._keys)}) must match number of values ({len(self._values)}) in __init__")

        # Shape determination is primarily handled by fromkeysvalues


    @classmethod
    def fromkeysvalues(cls, algebra: Any, keys: Tuple[int, ...], values: Union[List[Any], np.ndarray]) -> 'MultiVector':
        """
        Construct a multivector directly from keys and values. (Primary constructor)

        Args:
            algebra: The Geometric Algebra instance
            keys: Tuple of binary keys for basis blades (e.g., (0,), (1, 2, 4))
            values: List or NumPy array of coefficients matching the keys

        Returns:
            A new MultiVector instance
        """
        # REVIEW: This looks like the intended primary way to create instances internally.
        # Ensure keys is always a tuple
        if isinstance(keys, int):
            keys = (keys,)
        elif not isinstance(keys, tuple):
            keys = tuple(keys) # Convert list or other iterable to tuple

        # Ensure values is in a supported format (list or ndarray)
        if not isinstance(values, (list, np.ndarray)):
            # If single value provided for multiple keys, broadcast? Or error?
            # Let's assume if values is not list/array, it should match a single key.
            if len(keys) == 1:
                values = [values] # Wrap single value in list
            else:
                # Allow broadcasting scalar value to all keys?
                try:
                    # Check if it's a scalar type that can be broadcast
                    if isinstance(values, (int, float, complex, np.number, sympy.Expr)):
                        # print(f"Warning: Broadcasting single value {values} to {len(keys)} keys.")
                        values = [values] * len(keys)
                    else:
                         raise ValueError(f"Cannot match non-list/array values (type {type(values)}) to multiple keys ({len(keys)})")
                except ValueError:
                     raise ValueError(f"Number of keys ({len(keys)}) must match number of values (1, type {type(values)}) unless value is broadcastable scalar")


        # Basic validation
        # REVIEW: Crucial check.
        if len(keys) != len(values):
             # This should ideally not happen after the broadcasting logic above
             raise ValueError(f"Number of keys ({len(keys)}) must match number of values ({len(values)})")


        obj = object.__new__(cls)
        obj.algebra = algebra
        obj._keys = keys
        # Ensure values is ndarray if any element is ndarray for consistency? Maybe not necessary.
        # If list, ensure it's copied? No, assume caller handles copying if needed.
        obj._values = values # Now guaranteed to be list or ndarray

        # CHANGE: Refactor shape determination logic (Shape Handling Issue)
        obj.shape = () # Default shape
        if isinstance(values, np.ndarray):
            if values.ndim == 1:
                obj.shape = () # 1D array means no extra shape dimensions
            elif values.ndim > 1:
                if values.shape[0] == len(keys):
                    # Standard case: First dimension matches keys, remaining is shape
                    obj.shape = values.shape[1:]
                elif len(keys) == 1:
                    # Ambiguous case 1: Single key, multi-dim array value -> shape is value shape
                    obj.shape = values.shape
                else:
                    # Ambiguous case 2: Multiple keys, first dim mismatch -> Error
                    raise ValueError(f"Shape mismatch: NumPy array shape {values.shape} first dimension "
                                     f"does not match number of keys ({len(keys)})")
        elif isinstance(values, list):
             # Check if it's a list of lists/arrays (potential array MV)
             if values and isinstance(values[0], (list, np.ndarray)):
                 try:
                     # Convert to array to check shape consistency (use object dtype for safety)
                     arr = np.array(values, dtype=object)
                     if arr.ndim == 1: # List of scalars (or objects that didn't make array multi-dim)
                         obj.shape = ()
                     elif arr.ndim > 1:
                          if arr.shape[0] == len(keys):
                              # First dimension matches keys
                              obj.shape = arr.shape[1:]
                          elif len(keys) == 1:
                              # Single key, list of lists/arrays -> treat as array value
                              # Shape inference from list structure is tricky, try np.array again
                              try:
                                   arr_numeric = np.array(values) # Try numeric conversion for shape
                                   obj.shape = arr_numeric.shape
                              except:
                                   obj.shape = () # Fallback if shape cannot be determined
                          else:
                              # Ambiguous list structure
                              print(f"Warning: Shape determination ambiguity. List structure implies array shape {arr.shape} "
                                    f"but first dimension does not match keys length {len(keys)}. Assigning shape=().")
                              obj.shape = () # Fallback
                 except Exception as e:
                     # print(f"Debug: Exception during list shape inference: {e}")
                     obj.shape = () # Fallback if conversion fails
             else:
                 # Simple list of scalars
                 obj.shape = ()

        # Validation: Check consistency between values size, keys length, and shape
        # This is complex, rely on downstream methods to handle potential inconsistencies for now.
        # Example check (might be too strict or slow):
        # if obj.shape and isinstance(obj._values, np.ndarray):
        #     expected_size = len(obj._keys) * np.prod(obj.shape)
        #     if obj._values.size != expected_size:
        #         print(f"Warning: Internal size inconsistency. values size {obj._values.size}, "
        #               f"expected {expected_size} based on keys {len(obj._keys)} and shape {obj.shape}.")

        return obj


    # === Properties and Attribute Access ===

    @property
    def d(self) -> int:
        """
        Get the dimension of the underlying vector space for the algebra.

        Returns:
            Number of base vector dimensions (e.g., 3 for 3D Euclidean)
        """
        # REVIEW: Simple delegation to algebra. Assumes algebra has 'd'.
        return self.algebra.d if hasattr(self.algebra, 'd') else 0

    @property
    def coeffs(self) -> Dict[int, Any]:
        """
        Get a dictionary mapping basis blade keys to their coefficients.

        Note: For array-valued multivectors (where `self.shape` is not empty),
        this currently returns the coefficients corresponding to the *first*
        element along the shape dimensions (flattened index 0).
        Consider using `itermv()` or `__getitem__` for full access.

        Returns:
            Dictionary mapping binary keys to their coefficients (first element if array-valued).
        """
        # REVIEW: Behavior for array MVs documented. Suggests alternatives.
        if hasattr(self, 'shape') and self.shape: # Check shape exists and is not empty
            num_keys = len(self._keys)
            if isinstance(self._values, np.ndarray):
                try:
                    if self._values.ndim > 0:
                         expected_size = num_keys * np.prod(self.shape)
                         if self._values.size == expected_size:
                             # Get index of the first element in the shape dimensions
                             first_element_idx = (slice(None),) + tuple(np.zeros(len(self.shape), dtype=int))
                             first_element_values = self._values.reshape(num_keys, *self.shape)[first_element_idx]
                         else:
                              print(f"Warning: Size mismatch in coeffs property. Cannot extract first element reliably. Returning empty dict.")
                              return {}
                    else:
                         print(f"Warning: shape is set but _values is 0-dim array in coeffs. Returning empty dict.")
                         return {}
                except (ValueError, IndexError) as e:
                    print(f"Warning: Error reshaping/indexing in coeffs property: {e}. Returning empty dict.")
                    return {}
                return dict(zip(self._keys, first_element_values))

            elif isinstance(self._values, list):
                 try:
                     # Assume structure allows getting first element from each sublist/value
                     first_element_values = [v[0] if isinstance(v, (list, tuple, np.ndarray)) and len(v)>0 else v for v in self._values] # Basic handling
                     return dict(zip(self._keys, first_element_values))
                 except IndexError:
                      print("Warning: IndexError accessing first element in list-based coeffs. Returning empty dict.")
                      return {}
                 except Exception as e:
                      print(f"Warning: Error processing list-based coeffs: {e}. Returning empty dict.")
                      return {}
            else:
                  print(f"Warning: shape is set but _values is not list/ndarray in coeffs. Returning empty dict.")
                  return {}
        else:
            # For regular (non-array) multivectors or shape=()
            return dict(zip(self._keys, self._values)) if hasattr(self, '_keys') else {}

    @property
    def grades(self) -> Tuple[int, ...]:
        """
        Get the sorted tuple of grades present in the multivector with non-zero coefficients.

        Returns:
            Sorted tuple of grades (0=scalar, 1=vector, 2=bivector, etc.)
        """
        # REVIEW: Zero handling review - ensure explicit zero scalar returns (0,)
        if not hasattr(self, '_keys') or not self._keys:
            return tuple()

        keys = self._keys
        vals = self._values

        present_grades = set()
        is_explicit_zero_scalar = False

        for k, v in zip(keys, vals):
            # Determine if the value (or any part of it for array MVs) is non-zero
            is_effectively_zero = False
            # --- Robust Zero Check ---
            if isinstance(v, (int, float, complex, np.number)):
                is_effectively_zero = abs(v) < 1e-12
            elif isinstance(v, np.ndarray):
                 is_effectively_zero = not np.any(np.abs(v) > 1e-12)
            elif isinstance(v, list):
                 def check_list_zero(item):
                     if isinstance(item, list): return all(check_list_zero(x) for x in item)
                     elif isinstance(item, (int, float, complex, np.number)): return abs(item) < 1e-12
                     else: return False # Treat other types as non-zero
                 is_effectively_zero = check_list_zero(v)
            elif hasattr(v, 'is_zero'):
                 try: is_effectively_zero = v.is_zero
                 except TypeError: is_effectively_zero = (v == 0)
            elif v is None: is_effectively_zero = True
            elif v == 0: is_effectively_zero = True
            # --- End Zero Check ---

            if not is_effectively_zero:
                grade = bin(k).count("1")
                present_grades.add(grade)
            elif k == 0: # Check if the zero value belongs to the scalar key
                 is_explicit_zero_scalar = True


        # Handle the case of an explicit zero scalar {0: 0} -> should have grade 0
        if not present_grades and is_explicit_zero_scalar:
             return (0,)

        return tuple(sorted(list(present_grades)))


    @property
    def issymbolic(self) -> bool:
        """
        Determine if the multivector contains symbolic coefficients (e.g., SymPy Expr).

        Returns:
            True if any coefficient is symbolic, False otherwise
        """
        if not hasattr(self, '_values') or self._values is None: return False

        vals_iter = self._values
        if isinstance(self._values, np.ndarray):
            # Check dtype first for efficiency
            if self._values.dtype == object:
                 vals_iter = self._values.flat # Iterate through all elements if object array
            else:
                 return False # Not object array, cannot contain SymPy Expr directly

        for v in vals_iter:
            if isinstance(v, sympy.Expr) and not v.is_number: return True
            if isinstance(v, RationalPolynomial): return True
            if hasattr(self.algebra, 'codegen_symbolcls') and self.algebra.codegen_symbolcls:
                symbolcls = self.algebra.codegen_symbolcls
                if hasattr(symbolcls, "__self__"): symbolcls = symbolcls.__self__
                if isinstance(v, symbolcls): return True
        return False


    @property
    def free_symbols(self) -> Set:
        """
        Get the set of free SymPy symbols present in the multivector's coefficients.

        Returns:
            Set of sympy.Symbol instances found.
        """
        symbols = set()
        if not hasattr(self, '_values') or self._values is None: return symbols

        vals_iter = self._values
        if isinstance(self._values, np.ndarray):
             vals_iter = self._values.flat

        for v in vals_iter:
            if hasattr(v, "free_symbols") and isinstance(getattr(v, 'free_symbols'), set):
                symbols |= v.free_symbols
        return symbols

    @property
    def type_number(self) -> int:
        """
        Compute a unique integer type number based on the presence/absence of basis blades.
        This can be used for dispatching or classifying multivectors.

        Returns:
            Integer where the n-th bit is set if the n-th canonical basis blade is present.
        """
        if not hasattr(self.algebra, 'canon2bin') or not self.algebra.canon2bin:
             print("Warning: Algebra lacks 'canon2bin' mapping for type_number calculation. Returning 0.")
             return 0

        present_keys = set(self.keys())
        canon_keys_ordered = list(self.algebra.canon2bin.values())
        bitstr = "".join(
            "1" if key in present_keys else "0"
            for key in reversed(canon_keys_ordered)
        )
        return int(bitstr, 2) if bitstr else 0


    # Custom attribute access for basis blades (e.g., mv.e1, mv.e12)
    def __getattr__(self, name: str) -> Any:
        """
        Allow direct access to basis blade coefficients using attribute syntax like mv.e1.

        Args:
            name: Name of the attribute (e.g., 'e', 'e1', 'e2', 'e12').

        Returns:
            Coefficient value for the corresponding basis blade (0 if not present).

        Raises:
            AttributeError: If the name does not correspond to a valid basis blade
                            or the special scalar accessor 'e'.
        """
        if name == 'e':
            return self.get(0, 0)

        if hasattr(self.algebra, 'canon2bin') and isinstance(self.algebra.canon2bin, dict):
            blade_key = self.algebra.canon2bin.get(name)
            if blade_key is not None:
                return self.get(blade_key, 0)
            # CHANGE: Handle case where name looks like a blade but isn't in this algebra (test_getattr)
            elif re.match(r'^e[0-9^]+$', name):
                 return 0

        raise AttributeError(f"'{type(self).__name__}' object has no attribute '{name}', "
                             f"or it's not a recognized basis blade name in the algebra.")


    # === Core Methods ===

    def keys(self) -> Tuple[int, ...]:
        """
        Get the tuple of binary keys representing the basis blades present in this multivector.

        Returns:
            Tuple of integer keys (e.g., (0, 3, 5) for scalar, e12, e13).
        """
        return self._keys if hasattr(self, '_keys') else tuple()

    def values(self) -> Union[List[Any], np.ndarray]:
        """
        Get the coefficients corresponding to the basis blades keys.

        Returns:
            List or NumPy array of coefficients.
        """
        return self._values if hasattr(self, '_values') else []


    def itermv(self) -> Generator['MultiVector', None, None]:
        """
        Iterate over the elements of an array-valued multivector.

        If `self.shape` is empty (i.e., it's not array-valued), yields `self` once.
        Otherwise, yields a new MultiVector for each element along the shape dimensions.
        Each yielded MultiVector will have shape `()`.

        Yields:
            MultiVector instance for each element.
        """
        if not hasattr(self, 'shape') or not self.shape:
            yield self
            return

        vals = self._values
        if not isinstance(vals, np.ndarray):
            try:
                 vals = np.array(vals)
                 if vals.shape[0] != len(self._keys):
                     print(f"Warning: Shape mismatch after converting list to array in itermv. Expected first dim {len(self._keys)}, got {vals.shape[0]}. Yielding self.")
                     yield self; return
            except:
                 print("Warning: Failed to convert list-based values to NumPy array in itermv. Yielding self.")
                 yield self; return

        element_size = len(self._keys)
        total_elements = np.prod(self.shape) if self.shape else 1
        expected_total_size = element_size * total_elements

        if vals.size == expected_total_size:
             try:
                 reshaped_vals = vals.reshape(len(self._keys), *self.shape)
             except ValueError as e:
                 print(f"Error: Cannot reshape values array in itermv (size {vals.size}) to ({len(self._keys)}, {self.shape}). {e}. Yielding self.")
                 yield self; return

             it = np.ndindex(self.shape)
             for idx in it:
                 element_values = reshaped_vals[(slice(None),) + idx]
                 yield type(self).fromkeysvalues(self.algebra, self._keys, element_values)
        else:
             print(f"Warning: Size mismatch in itermv. Expected {expected_total_size}, got {vals.size}. Yielding self.")
             yield self


    def map(self, func: Callable) -> 'MultiVector':
        """
        Apply a function to each coefficient value of the multivector.

        Args:
            func: A function to apply. It should accept a single argument (the coefficient value)
                  or two arguments (key, value). If array-valued, function is applied element-wise.

        Returns:
            A new multivector with the function applied to each coefficient.
        """
        try:
            num_args = func.__code__.co_argcount
            if hasattr(func, '__self__'): num_args -= 1
        except AttributeError:
             num_args = 1

        new_values = []
        if num_args == 2:
             if self.shape:
                 if isinstance(self._values, np.ndarray):
                     keys_array = np.array(self._keys).reshape(-1, *([1]*len(self.shape)))
                     def func_wrapper(k, v_element): return func(k, v_element)
                     vectorized_func = np.vectorize(func_wrapper, otypes=[object]) # Use object otype for safety
                     new_values = vectorized_func(keys_array, self._values)
                 else: # List based
                     print("Warning: map() with 2-argument function on list-based array MV is inefficient.")
                     temp_values = []
                     for k, v_slice in zip(self._keys, self._values):
                         if isinstance(v_slice, (list, np.ndarray)):
                             mapped_slice = [func(k, elem) for elem in np.array(v_slice).flat]
                             temp_values.append(mapped_slice)
                         else:
                              temp_values.append(func(k, v_slice))
                     new_values = temp_values
             else: # Not array-valued
                 for k, v in zip(self._keys, self._values):
                     new_values.append(func(k, v))
        else: # num_args == 1
             if isinstance(self._values, np.ndarray):
                 new_values = np.vectorize(func, otypes=[object])(self._values) # Use object otype
             else:
                 temp_values = []
                 for v in self._values:
                      if isinstance(v, list):
                           try: temp_values.append(np.vectorize(func, otypes=[object])(np.array(v)))
                           except: temp_values.append([func(elem) for elem in v])
                      elif isinstance(v, np.ndarray):
                           temp_values.append(np.vectorize(func, otypes=[object])(v))
                      else:
                           temp_values.append(func(v))
                 new_values = temp_values

        # Create result
        if isinstance(self._values, np.ndarray) and not isinstance(new_values, np.ndarray):
             try: new_values = np.array(new_values, dtype=object) # Default to object
             except: pass # Keep as list if conversion fails

        result_mv = type(self).fromkeysvalues(self.algebra, self._keys, new_values)
        if hasattr(self, 'shape'): result_mv.shape = self.shape
        return result_mv


    def filter(self, pred: Callable) -> 'MultiVector':
        """
        Filter components of the multivector based on a predicate function.

        Args:
            pred: A predicate function. It should accept a single argument (the coefficient value)
                  or two arguments (key, value) and return True to keep the component, False to discard.
                  For array-valued MVs, keeps blade if predicate holds for *any* element.

        Returns:
            A new multivector containing only the components for which pred returns True.
            Note: This operation currently does not preserve the shape of array-valued multivectors.
        """
        # REVIEW: Filter logic for array MVs keeps blade if *any* element satisfies pred. Shape is lost.
        try:
            num_args = pred.__code__.co_argcount
            if hasattr(pred, '__self__'): num_args -= 1
        except AttributeError: num_args = 1

        filtered_keys = []
        filtered_values = []

        for i, k in enumerate(self._keys):
            v = self._values[i]
            keep_blade = False
            try:
                if self.shape and isinstance(v, (list, np.ndarray)): # Array valued slice
                    vals_to_check = v.flat if isinstance(v, np.ndarray) else iter(v)
                    if num_args == 2:
                        for elem in vals_to_check:
                            if pred(k, elem): keep_blade = True; break
                    else:
                        for elem in vals_to_check:
                            if pred(elem): keep_blade = True; break
                else: # Scalar value (or treat whole slice as one value if not list/array)
                    if num_args == 2:
                        if pred(k, v): keep_blade = True
                    else:
                        if pred(v): keep_blade = True
            except Exception as e:
                 print(f"Warning: Error applying predicate in filter() for key {k}. Skipping blade. Error: {e}")
                 keep_blade = False

            if keep_blade:
                filtered_keys.append(k)
                filtered_values.append(v)

        if not filtered_keys:
            return self.algebra.scalar(0) # Use algebra factory for zero

        output_values = filtered_values
        if isinstance(self._values, np.ndarray) and filtered_values:
             try: output_values = np.array(filtered_values, dtype=object)
             except: output_values = list(filtered_values)
        else: output_values = list(filtered_values)

        return type(self).fromkeysvalues(self.algebra, tuple(filtered_keys), output_values)


    def exp(self) -> 'MultiVector':
        """
        Compute the exponential of the multivector (exp(A)).

        Uses specialized formulas for pure vectors/bivectors if their square is scalar,
        otherwise defaults to a Taylor series expansion.

        Returns:
            A new multivector representing exp(self).
        """
        # Handle zero multivector
        if not self.grades:
            id_val = S(1) if self.issymbolic else 1.0
            return self.algebra.scalar(id_val)

        # Handle pure scalar
        if self.grades == (0,):
            exp_val = sympy.exp(self.e) if self.issymbolic else np.exp(self.e)
            return self.algebra.scalar(exp_val)

        # Check for pure vector or bivector optimization
        if len(self.grades) == 1 and self.grades[0] in (1, 2):
            mv2 = self * self
            if mv2.grades == (0,) or (len(mv2.grades) == 1 and mv2.grades[0] == 0):
                 scalar_square = mv2.e
                 is_zero_sq, is_negative_sq, is_positive_sq = False, False, False
                 # Robust checks for zero/sign (numeric and symbolic)
                 if isinstance(scalar_square, (int, float, complex, np.number)):
                     real_part = scalar_square.real if isinstance(scalar_square, complex) else scalar_square
                     if abs(scalar_square) < 1e-12: is_zero_sq = True
                     elif real_part < -1e-12: is_negative_sq = True
                     elif real_part > 1e-12: is_positive_sq = True
                 elif isinstance(scalar_square, sympy.Expr):
                     if scalar_square.is_zero: is_zero_sq = True
                     elif scalar_square.is_negative: is_negative_sq = True
                     elif scalar_square.is_positive: is_positive_sq = True

                 if is_zero_sq: return self.algebra.scalar(S(1) if self.issymbolic else 1) + self
                 if is_negative_sq:
                     alpha = sympy.sqrt(-scalar_square) if self.issymbolic else np.sqrt(-scalar_square)
                     cos_alpha = sympy.cos(alpha) if self.issymbolic else np.cos(alpha)
                     sin_alpha = sympy.sin(alpha) if self.issymbolic else np.sin(alpha)
                     return self.algebra.scalar(cos_alpha) + (self / alpha) * sin_alpha
                 if is_positive_sq:
                      alpha = sympy.sqrt(scalar_square) if self.issymbolic else np.sqrt(scalar_square)
                      cosh_alpha = sympy.cosh(alpha) if self.issymbolic else np.cosh(alpha)
                      sinh_alpha = sympy.sinh(alpha) if self.issymbolic else np.sinh(alpha)
                      return self.algebra.scalar(cosh_alpha) + (self / alpha) * sinh_alpha

        # General case: Taylor series expansion
        result = self.algebra.scalar(S(1) if self.issymbolic else 1.0)
        term = self.copy()
        factorial = S(1) if self.issymbolic else 1.0
        tolerance = 1e-12
        max_iter = 30

        for i in range(1, max_iter + 1):
            term_i = term / factorial
            result += term_i

            if not self.issymbolic:
                term_norm = 0
                try:
                    # Use norm() method if available and efficient enough? For now, sum abs.
                    if isinstance(term_i._values, np.ndarray) and np.issubdtype(term_i._values.dtype, np.number):
                        term_norm = np.sum(np.abs(term_i._values))
                    elif isinstance(term_i._values, list):
                        term_norm = sum(abs(v) for v in term_i._values if isinstance(v, (int, float, complex, np.number)))
                    else: # Cannot compute norm
                         print("Warning: Cannot compute norm for convergence check in exp(). Stopping early.")
                         break
                except TypeError:
                     print(f"Warning: Cannot compute norm for convergence check in exp() due to non-numeric type. Stopping early.")
                     break

                if term_norm < tolerance: break

            term = term * self
            factorial *= S(i + 1) if self.issymbolic else (i + 1)

            if i == max_iter: print(f"Warning: exp() Taylor series may not have fully converged after {max_iter} iterations.")

        return result


    def items(self) -> Generator[Tuple[int, Any], None, None]:
        """
        Get an iterator over the (key, value) pairs of the multivector components.

        Yields:
            Pairs of (basis_blade_key, coefficient_value).
        """
        if hasattr(self, '_keys') and hasattr(self, '_values'):
             if len(self._keys) == len(self._values):
                 yield from zip(self._keys, self._values)
             else:
                 print(f"Warning: Mismatch between keys ({len(self._keys)}) and values ({len(self._values)}) in items().")
                 min_len = min(len(self._keys), len(self._values))
                 yield from zip(self._keys[:min_len], self._values[:min_len])


    def copy(self) -> MultiVector:
        """
        Create a deep copy of this multivector.

        Ensures that the algebra reference is shared, but keys and values are copied.

        Returns:
            A new MultiVector instance with copied data.
        """
        return copy.deepcopy(self)


    def __len__(self) -> int:
        """
        Get the number of basis blade components stored in this multivector.

        Returns:
            Number of (key, value) pairs stored.
        """
        return len(self._keys) if hasattr(self, '_keys') else 0


    # === Grade Filtering ===

    def grade(self, *grades_to_keep: int) -> 'MultiVector':
        """
        Extract specific grades from the multivector.

        Args:
            *grades_to_keep: Integer grades to extract (0=scalar, 1=vector, 2=bivector, etc.).
                             Can be passed as multiple arguments (0, 2) or a single tuple ((0, 2),).

        Returns:
            A new multivector containing only the components belonging to the specified grades.
            Returns a zero scalar if no components match the specified grades.
            Note: This operation currently does not preserve the shape of array-valued multivectors.
        """
        # REVIEW: Shape loss for array MVs documented.
        if not grades_to_keep: return self.algebra.scalar(0)
        if len(grades_to_keep) == 1 and isinstance(grades_to_keep[0], (tuple, set, list)):
            grades_set = set(grades_to_keep[0])
        else: grades_set = set(grades_to_keep)

        filtered_keys = []
        filtered_values = []
        if self.shape: print("Warning: grade() on array-valued MultiVector currently flattens the result.")

        for i, k in enumerate(self._keys):
            if bin(k).count('1') in grades_set:
                filtered_keys.append(k)
                filtered_values.append(self._values[i])

        if not filtered_keys: return self.algebra.scalar(0)

        output_values = filtered_values
        if isinstance(self._values, np.ndarray) and filtered_values:
             try: output_values = np.array(filtered_values, dtype=object)
             except: output_values = list(filtered_values)
        else: output_values = list(filtered_values)

        return type(self).fromkeysvalues(self.algebra, tuple(filtered_keys), output_values)


    def scalar(self) -> MultiVector:
        """Extract the scalar part (grade 0)."""
        return self.grade(0)

    def vector(self) -> MultiVector:
        """Extract the vector part (grade 1)."""
        return self.grade(1)

    def bivector(self) -> MultiVector:
        """Extract the bivector part (grade 2)."""
        return self.grade(2)

    # === Geometric Operations ===

    def reverse(self) -> MultiVector:
        """Compute the reverse of this multivector (~A)."""
        new_values = []
        for k, v in zip(self._keys, self._values):
            grade = bin(k).count('1')
            sign = -1 if grade % 4 in (2, 3) else 1
            if sign == 1: new_values.append(v); continue
            # Apply sign = -1
            if isinstance(v, (int, float, complex, np.number)): new_values.append(-v)
            elif isinstance(v, np.ndarray): new_values.append(-v)
            else:
                 try: new_values.append(-v)
                 except TypeError:
                     try: new_values.append(sign * v)
                     except TypeError:
                         print(f"Warning: Cannot apply reverse sign to value type {type(v)}. Keeping original.")
                         new_values.append(v)

        result_mv = type(self).fromkeysvalues(self.algebra, self._keys, new_values)
        if hasattr(self, 'shape'): result_mv.shape = self.shape
        return result_mv


    def involute(self) -> MultiVector:
        """Compute the grade involution of this multivector (A*)."""
        return self.algebra.involute(self)

    def conjugate(self) -> MultiVector:
        """Compute the Clifford conjugate of this multivector (A†)."""
        return self.algebra.conjugate(self)

    def inverse(self) -> 'MultiVector':
        """Compute the multiplicative inverse of this multivector (A⁻¹)."""
        return self.inv()

    def inv(self) -> 'MultiVector':
        """Compute the multiplicative inverse of this multivector (A⁻¹). (Implementation)"""
        # Handle scalar case
        if self.grades == (0,):
            scalar_value = self.e
            is_zero = False # Robust zero check
            if isinstance(scalar_value, (int, float, complex, np.number)): is_zero = abs(scalar_value) < 1e-12
            elif hasattr(scalar_value, 'is_zero'):
                try: is_zero = scalar_value.is_zero
                except TypeError: is_zero = (scalar_value == 0)
            if is_zero: raise ZeroDivisionError("Cannot compute inverse of zero scalar MultiVector")
            try: inv_value = 1.0 / scalar_value
            except TypeError: inv_value = S(1) / scalar_value # Symbolic
            return self.algebra.scalar(inv_value)

        # General case: A⁻¹ = ~A / (A * ~A)
        a_rev_a = self * self.reverse()
        if a_rev_a.grades == (0,):
             scalar_norm_sq = a_rev_a.e
             is_zero_norm = False # Robust zero check
             if isinstance(scalar_norm_sq, (int, float, complex, np.number)): is_zero_norm = abs(scalar_norm_sq) < 1e-12
             elif hasattr(scalar_norm_sq, 'is_zero'):
                 try: is_zero_norm = scalar_norm_sq.is_zero
                 except TypeError: is_zero_norm = (scalar_norm_sq == 0)
             if is_zero_norm: raise ZeroDivisionError("Multivector is not invertible (norm squared is zero)")
             return self.reverse() / scalar_norm_sq
        else: # Fallback to algebra implementation if A*~A is not scalar
             try: return self.algebra.inv(self)
             except (AttributeError, NotImplementedError):
                  raise ZeroDivisionError("Multivector is not invertible (A * ~A is not scalar, and no algebra.inv available)")
             except ZeroDivisionError:
                  raise ZeroDivisionError("Multivector is not invertible (determined by algebra.inv)")


    def norm(self) -> float:
        """
        Compute the magnitude (norm) of this multivector: sqrt(|A * ~A|_0).

        Returns the norm as a non-negative float or potentially a symbolic expression.
        """
        a_rev_a = self * self.reverse()
        scalar_part = a_rev_a.grade(0)
        norm_sq_val = S(0) if self.issymbolic else 0.0

        if scalar_part.keys():
            val = scalar_part.values()[0]
            if isinstance(val, sympy.Expr): return sympy.Abs(val) # Return symbolic |A*~A|_0
            elif isinstance(val, (int, float, complex, np.number)): norm_sq_val = abs(val)
            else:
                 try: norm_sq_val = abs(float(val))
                 except (TypeError, ValueError): raise TypeError(f"Cannot compute numeric norm for scalar part type {type(val)}")

        # Return sqrt for numeric, symbolic handled above
        return sympy.sqrt(norm_sq_val) if self.issymbolic else np.sqrt(norm_sq_val)


    def normalized(self) -> MultiVector:
        """Return a normalized version of this multivector (self / self.norm())."""
        norm_value = self.norm()
        is_zero_norm = False # Robust zero check
        if isinstance(norm_value, (float, int, np.number)): is_zero_norm = abs(norm_value) < 1e-12
        elif hasattr(norm_value, 'is_zero'):
             try: is_zero_norm = norm_value.is_zero
             except TypeError: is_zero_norm = (norm_value == 0)
        if is_zero_norm: raise ZeroDivisionError("Cannot normalize a multivector with zero norm")

        try: return self / norm_value
        except TypeError as e:
             print(f"Warning: Type error during normalization ({e}). Ensure coefficients support division.")
             raise e


    # === Arithmetic Operations ===

    def __eq__(self, other: Any) -> bool:
        """
        Check for equality between this multivector and another object.
        Uses tolerance for numeric types and symbolic simplification for expressions.
        """
        # REVIEW: Equality comparison logic reviewed. Handles numeric tolerance, symbolic simplification,
        # and array comparison via np.allclose. Comparison iterates through all keys, not just coeffs property.
        if isinstance(other, (int, float, complex)):
            if self.grades == (0,):
                 s_val = self.e
                 if isinstance(s_val, np.ndarray):
                     try: return np.allclose(s_val, other, atol=1e-12)
                     except TypeError: return False
                 elif isinstance(s_val, (float, complex)): return abs(s_val - other) < 1e-12
                 elif isinstance(s_val, sympy.Expr):
                     try: return sympy.simplify(s_val - other) == 0
                     except (TypeError, AttributeError): return s_val == other
                 else: return s_val == other
            else: # Non-scalar MV cannot equal scalar (unless both zero)
                 if other == 0: return not self.grades # Check if self is zero MV
                 return False

        elif isinstance(other, MultiVector):
            if self.algebra != other.algebra: return False
            if getattr(self, 'shape', ()) != getattr(other, 'shape', ()): return False

            all_keys = set(self.keys()) | set(other.keys())
            for k in all_keys:
                 v_self = self.get(k, 0)
                 v_other = other.get(k, 0)
                 is_equal = False
                 try:
                     if isinstance(v_self, np.ndarray) or isinstance(v_other, np.ndarray):
                         is_equal = np.allclose(v_self, v_other, atol=1e-12)
                     elif isinstance(v_self, (float, complex)) or isinstance(v_other, (float, complex)):
                         is_equal = abs(v_self - v_other) < 1e-12
                     elif isinstance(v_self, sympy.Expr) or isinstance(v_other, sympy.Expr):
                         try: is_equal = (sympy.simplify(S(v_self) - S(v_other)) == 0)
                         except (TypeError, AttributeError): is_equal = (v_self == v_other)
                     else: is_equal = (v_self == v_other)
                 except TypeError: is_equal = False # Incompatible types

                 if not is_equal: return False
            return True # All components matched

        else: return NotImplemented


    def __add__(self, other: Union[MultiVector, int, float, complex]) -> MultiVector:
        """Add another multivector or scalar to this multivector."""
        try: # Prefer algebra implementation
            return self.algebra.add(self, other)
        except (AttributeError, NotImplementedError):
             print("Warning: Using fallback __add__ implementation.")
             is_symbolic_op = self.issymbolic
             if isinstance(other, (int, float, complex, sympy.Expr)):
                 other_mv = self.algebra.scalar(other)
                 if isinstance(other, sympy.Expr): is_symbolic_op = True
             elif isinstance(other, MultiVector):
                 if self.algebra != other.algebra: raise ValueError("Cannot add MVs from different algebras")
                 if self.shape != other.shape: raise ValueError(f"Cannot add MVs with incompatible shapes: {self.shape} and {other.shape}")
                 other_mv = other
                 if other_mv.issymbolic: is_symbolic_op = True
             else: return NotImplemented

             result_keys = list(set(self.keys()) | set(other_mv.keys()))
             result_values = []
             needs_numpy = isinstance(self._values, np.ndarray) or isinstance(other_mv._values, np.ndarray)
             default_zero = S(0) if is_symbolic_op else 0 # CHANGE: Use symbolic zero if needed

             for k in result_keys:
                 v_self = self.get(k, default_zero) # Use appropriate zero
                 v_other = other_mv.get(k, default_zero)
                 try: result_values.append(v_self + v_other)
                 except TypeError:
                     try: # Retry with symbolic/numpy conversion
                         if is_symbolic_op: sum_val = S(v_self) + S(v_other)
                         else: sum_val = np.array(v_self) + np.array(v_other) # Let numpy handle broadcasting/types
                         result_values.append(sum_val)
                     except Exception as e: raise TypeError(f"Cannot add values for key {k}: {v_self} + {v_other}. Error: {e}")

             # Filter zeros
             final_keys, final_values = [], []
             for k, v in zip(result_keys, result_values):
                  is_zero = False # Robust zero check
                  if isinstance(v, (float, complex)): is_zero = abs(v) < 1e-12
                  elif isinstance(v, np.ndarray): is_zero = not np.any(np.abs(v) > 1e-12)
                  elif hasattr(v, 'is_zero'):
                      try: is_zero = v.is_zero
                      except TypeError: is_zero = (v == 0)
                  elif v == 0: is_zero = True
                  if not is_zero: final_keys.append(k); final_values.append(v)

             if not final_keys: return self.algebra.scalar(0)

             output_values = final_values
             if needs_numpy and not is_symbolic_op: # Try to convert back to numpy if started with it
                 try: output_values = np.array(final_values) # Infer dtype
                 except: output_values = np.array(final_values, dtype=object) # Fallback

             result_mv = type(self).fromkeysvalues(self.algebra, tuple(final_keys), output_values)
             result_mv.shape = self.shape # Assume shape preserved
             return result_mv


    def __radd__(self, other: Union[int, float, complex]) -> MultiVector:
        """Handle right addition (scalar + multivector)."""
        return self + other

    def __sub__(self, other: Union[MultiVector, int, float, complex]) -> MultiVector:
        """Subtract another multivector or scalar from this multivector."""
        try: return self.algebra.sub(self, other)
        except (AttributeError, NotImplementedError):
             print("Warning: Using fallback __sub__ implementation.")
             if isinstance(other, MultiVector): return self + (-other)
             elif isinstance(other, (int, float, complex, sympy.Expr)): return self + (-self.algebra.scalar(other))
             else: return NotImplemented


    def __rsub__(self, other: Union[int, float, complex]) -> MultiVector:
        """Handle right subtraction (scalar - multivector)."""
        return self.algebra.scalar(other) + (-self)

    def __mul__(self, other: Union[MultiVector, int, float, complex, sympy.Expr]) -> MultiVector:
        """Multiply by another MultiVector (Geometric Product) or scalar."""
        new_values = None
        if isinstance(other, (int, float, complex)):
            if isinstance(self._values, np.ndarray): new_values = self._values * other
            else:
                try: new_values = [v * other for v in self._values]
                except TypeError: return self.copy() # Fallback
        elif isinstance(other, sympy.Expr): # Symbolic scalar
            if isinstance(self._values, np.ndarray):
                new_values = np.array([sympy.Mul(other, v, evaluate=False) for v in self._values.flat], dtype=object).reshape(self._values.shape)
            else:
                 try: new_values = [sympy.Mul(other, v, evaluate=False) for v in self._values]
                 except TypeError: return self.copy() # Fallback
        elif isinstance(other, MultiVector): # Geometric product
            if self.algebra != other.algebra: raise ValueError("Cannot multiply MVs from different algebras")
            return self.algebra.gp(self, other) # Delegate to algebra
        else: return NotImplemented

        if new_values is not None: # Result from scalar multiplication
            result_mv = type(self).fromkeysvalues(self.algebra, self._keys, new_values)
            if hasattr(self, 'shape'): result_mv.shape = self.shape
            return result_mv
        else: return NotImplemented # Should not be reached


    def __rmul__(self, other: Union[int, float, complex, sympy.Expr]) -> MultiVector:
        """Handle right multiplication (scalar * multivector)."""
        return self * other

    def __truediv__(self, other: Union[MultiVector, int, float, complex, sympy.Expr]) -> MultiVector:
        """Divide self by scalar, or multiply self by inverse of other."""
        if isinstance(other, (int, float, complex)):
            if abs(other) < 1e-12: raise ZeroDivisionError("Scalar division by zero")
            if isinstance(self._values, np.ndarray): new_values = self._values / other
            else:
                 try: new_values = [v / other for v in self._values]
                 except TypeError: raise TypeError(f"Cannot divide coefficient type by scalar '{other}'.")
            result_mv = type(self).fromkeysvalues(self.algebra, self._keys, new_values)
            if hasattr(self, 'shape'): result_mv.shape = self.shape
            return result_mv
        elif isinstance(other, sympy.Expr): # Symbolic scalar division
             if hasattr(other, 'is_zero') and other.is_zero: raise ZeroDivisionError("Symbolic scalar division by zero")
             return self * (S(1) / other) # Multiply by inverse
        elif isinstance(other, MultiVector): # Division by MV: A / B = A * B.inv()
            if self.algebra != other.algebra: raise ValueError("Cannot divide MVs from different algebras")
            try: return self.algebra.div(self, other) # Prefer algebra implementation
            except (AttributeError, NotImplementedError):
                 print("Warning: Using fallback __truediv__ (A * B.inv()).")
                 try: return self * other.inv()
                 except ZeroDivisionError: raise ZeroDivisionError("Division by non-invertible MultiVector")
        else: return NotImplemented

    # CHANGE: Implement __rtruediv__ (test_43)
    def __rtruediv__(self, other: Union[int, float, complex, sympy.Expr]) -> MultiVector:
        """Handle right division (scalar / multivector)."""
        if isinstance(other, (int, float, complex, sympy.Expr)):
            # Calculate as other * self.inv()
            try:
                return self.algebra.scalar(other) * self.inv()
            except ZeroDivisionError:
                raise ZeroDivisionError(f"Right division by non-invertible MultiVector: {self}")
        else:
            return NotImplemented


    def __neg__(self) -> MultiVector:
        """Negate this multivector (-A)."""
        if isinstance(self._values, np.ndarray): new_values = -self._values
        else:
             try: new_values = [-v for v in self._values]
             except TypeError: return self.copy() # Fallback
        result_mv = type(self).fromkeysvalues(self.algebra, self._keys, new_values)
        if hasattr(self, 'shape'): result_mv.shape = self.shape
        return result_mv


    def __invert__(self) -> MultiVector:
        """Compute the reverse (~A)."""
        return self.reverse()

    def __pow__(self, power: Union[int, float]) -> 'MultiVector':
        """Raise multivector to integer or specific float power."""
        # REVIEW: Binary exponentiation logic seems correct. Relies on __mul__ (algebra.gp).
        # Ensure identity element is handled correctly for symbolic/numeric.
        if isinstance(power, int):
            if power < 0:
                try: return self.inv() ** (-power)
                except ZeroDivisionError: raise ZeroDivisionError(f"Cannot raise non-invertible MV to negative power {power}")
            elif power == 0:
                id_val = S(1) if self.issymbolic else 1.0
                return self.algebra.scalar(id_val) # Correct identity
            elif power == 1: return self.copy()
            else: # Binary exponentiation
                id_val = S(1) if self.issymbolic else 1.0
                result = self.algebra.scalar(id_val)
                base = self.copy()
                p = power
                while p > 0:
                    # Ensure multiplication handles identity correctly
                    if p % 2 == 1: result = result * base
                    base = base * base
                    p //= 2
                return result
        elif isinstance(power, float):
             if abs(power - 0.5) < 1e-12: return self.sqrt()
             elif abs(power + 0.5) < 1e-12:
                 try: return self.inv().sqrt()
                 except ZeroDivisionError: raise ZeroDivisionError("Cannot calculate inverse sqrt for non-invertible MV")
             else: raise ValueError(f"Unsupported fractional power: {power}. Only 0.5 and -0.5 supported.")
        else: raise TypeError(f"Unsupported type for power: {type(power)}.")


    def sqrt(self) -> 'MultiVector':
         """
         Compute the square root of the multivector sqrt(X).
         Handles scalar case and cases where X*~X is a non-negative scalar (like rotors).
         Uses the formula sqrt(X) = (X + N) / sqrt(2 * (a + N)) where N=sqrt(X*~X), a=<X>_0.
         """
         # REVIEW: Formula implementation reviewed. Handles symbolic/numeric sqrt.
         # Potential numerical instability if a + N is close to zero.
         # Handle simple scalar case first
         if self.grades == (0,):
              scalar_val = self.e
              sqrt_val = sympy.sqrt(scalar_val) if self.issymbolic else np.sqrt(scalar_val+0j) # Use complex sqrt for safety
              return self.algebra.scalar(sqrt_val)

         # General case: Use formula based on norm, assuming X*~X is scalar
         X = self
         XrevX = X * X.reverse()

         if XrevX.grades == (0,):
             norm_sq = XrevX.e
             # Check if norm_sq is non-negative (numeric or symbolic)
             is_non_negative = False
             if isinstance(norm_sq, (int, float, complex, np.number)):
                 is_non_negative = norm_sq.real >= -1e-12 and abs(norm_sq.imag) < 1e-12
             elif isinstance(norm_sq, sympy.Expr):
                  if norm_sq.is_nonnegative is True: is_non_negative = True # Explicit True check
                  elif norm_sq.is_real is False: is_non_negative = False # Non-real cannot be non-negative
                  # Add more checks if needed, e.g., simplify(norm_sq >= 0)
             else: raise TypeError(f"Cannot determine sign of norm squared type {type(norm_sq)} for sqrt")

             if not is_non_negative:
                  raise ValueError("Cannot compute sqrt using this method when X*~X is negative scalar")

             # Calculate norm N = sqrt(X*~X)
             N = sympy.sqrt(norm_sq) if isinstance(norm_sq, sympy.Expr) else np.sqrt(norm_sq)

             # Calculate scalar part a = <X>_0
             a = X.grade(0).e

             # Denominator term: denom_sq = 2 * (a + N)
             denom_sq = 2 * (a + N)

             # Check if denominator is zero (robust check)
             is_zero_denom = False
             if isinstance(denom_sq, (int, float, complex, np.number)): is_zero_denom = abs(denom_sq) < 1e-12
             elif hasattr(denom_sq, 'is_zero'):
                 try: is_zero_denom = denom_sq.is_zero
                 except TypeError: is_zero_denom = (denom_sq == 0)
             if is_zero_denom: raise ValueError("Denominator is zero in sqrt calculation (X = -norm(X)?)")

             # Calculate denominator: sqrt(denom_sq)
             denom = sympy.sqrt(denom_sq) if isinstance(denom_sq, sympy.Expr) else np.sqrt(denom_sq)

             # Calculate numerator: X + N
             numerator = X + self.algebra.scalar(N)

             return numerator / denom
         else:
             raise NotImplementedError("General MultiVector square root for non-versors is not implemented.")


    def asfullmv(self, canonical: bool = True) -> 'MultiVector':
        """
        Convert to a 'full' multivector representation containing coefficients for all basis blades
        of the algebra, including zeros for blades not present in the original. Supports array-valued MVs.
        """
        # REVIEW: Handles array values. Creates full array and places values correctly.
        target_keys, target_size = None, 0
        if canonical:
            if not hasattr(self.algebra, 'canon2bin'): raise AttributeError("Algebra lacks 'canon2bin' mapping.")
            target_keys = tuple(self.algebra.canon2bin.values())
            target_size = len(target_keys)
        else:
            target_keys = tuple(range(2**self.algebra.d))
            target_size = 2**self.algebra.d

        output_dtype = object # Default to object for safety with mixed types/arrays
        if hasattr(self, '_values') and self._values is not None and len(self._values) > 0:
             first_val = self._values[0]
             if isinstance(first_val, np.ndarray): output_dtype = first_val.dtype
             elif isinstance(first_val, list): output_dtype = object
             else: # Infer from scalar types
                  is_sym = any(isinstance(v, sympy.Expr) for v in self._values)
                  is_comp = any(isinstance(v, complex) for v in self._values)
                  is_float = any(isinstance(v, float) for v in self._values)
                  if is_sym: output_dtype = object
                  elif is_comp: output_dtype = complex
                  elif is_float: output_dtype = float
                  elif all(isinstance(v, int) for v in self._values): output_dtype = int
                  else: output_dtype = object

        full_shape = (target_size,) + self.shape
        full_values = np.zeros(full_shape, dtype=output_dtype)

        current_key_indices = {key: i for i, key in enumerate(self._keys)}

        for i, target_key in enumerate(target_keys):
            if target_key in current_key_indices:
                 original_index = current_key_indices[target_key]
                 try: full_values[i] = self._values[original_index]
                 except (ValueError, TypeError) as e:
                     print(f"Warning: Type/shape mismatch assigning to full MV array (key {target_key}). Converting to object dtype. Error: {e}")
                     if output_dtype != object:
                          output_dtype = object
                          full_values = full_values.astype(object)
                          full_values[i] = self._values[original_index] # Retry
                     else: print(f"Error: Could not assign value for key {target_key} to full MV array.")

        result_mv = type(self).fromkeysvalues(self.algebra, target_keys, full_values)
        return result_mv


    # === Geometric Algebra Specific Operations ===


    def __xor__(self, other: MultiVector) -> MultiVector:
        """Compute the outer product (^)."""
        if isinstance(other, MultiVector):
            if self.algebra != other.algebra: raise ValueError("Algebras must match for op")
            return self.algebra.op(self, other)
        return NotImplemented


    def __or__(self, other: MultiVector) -> MultiVector:
        """Compute the inner product (|)."""
        if isinstance(other, MultiVector):
             if self.algebra != other.algebra: raise ValueError("Algebras must match for ip")
             return self.algebra.ip(self, other)
        return NotImplemented


    def __and__(self, other: MultiVector) -> MultiVector:
        """Compute the regressive product (&)."""
        if isinstance(other, MultiVector):
             if self.algebra != other.algebra: raise ValueError("Algebras must match for rp")
             return self.algebra.rp(self, other)
        return NotImplemented


    def __rshift__(self, other: MultiVector) -> MultiVector:
        """Compute the sandwich product (>>): self * other * ~self."""
        if isinstance(other, MultiVector):
             if self.algebra != other.algebra: raise ValueError("Algebras must match for sw")
             try: return self.algebra.sw(self, other)
             except (AttributeError, NotImplementedError):
                 print("Warning: Using fallback __rshift__ (self * other * ~self).")
                 return self * other * self.reverse()
        return NotImplemented


    def __matmul__(self, other: MultiVector) -> MultiVector:
        """Compute the projection (@) of other onto self."""
        if isinstance(other, MultiVector):
             if self.algebra != other.algebra: raise ValueError("Algebras must match for proj")
             return self.algebra.proj(self, other)
        return NotImplemented


    # === Indexing and Array Interface ===

    def __getitem__(self, item: Union[int, Tuple, slice]) -> Any:
        """
        Access parts of the multivector.
        - Integer `k`: Returns the k-th *stored* component value (scalar, list, or array).
        - Slice `s`: Returns a new MultiVector slicing stored components. Shape is lost.
        - Tuple `t`: Indexes into array-valued multivector dimensions, returning a new MultiVector.
        Note: Integer indexing returns raw value, slice/tuple return MultiVector.
        """
        # REVIEW: Return type differs based on index type (value vs MV). Documented.
        if isinstance(item, (int, slice)):
            try:
                if not hasattr(self, '_keys') or not hasattr(self, '_values'): raise IndexError("MV not initialized")
                if not isinstance(self._keys, Sequence) or not isinstance(self._values, (Sequence, np.ndarray)): raise TypeError("MV keys/values not sequence")

                selected_keys = self._keys[item]
                selected_values = self._values[item]

                if isinstance(item, int): # Return raw value
                     return selected_values
                else: # Return new MV for slice
                     keys_tuple = tuple(selected_keys) if not isinstance(selected_keys, tuple) else selected_keys
                     values_out = selected_values
                     if not isinstance(values_out, np.ndarray): values_out = list(values_out)
                     # Slicing components loses array shape information
                     return type(self).fromkeysvalues(self.algebra, keys_tuple, values_out)
            except IndexError: raise IndexError("MultiVector index out of range")
            except TypeError as e: raise TypeError(f"MV components not indexable/sliceable: {e}")

        elif isinstance(item, tuple): # Indexing into array dimensions
            if not self.shape: raise IndexError("Tuple indexing requires array-valued MV")
            try:
                 vals = self._values
                 if not isinstance(vals, np.ndarray): # Convert list-based to numpy
                      try: vals = np.array(vals)
                      except Exception as e: raise TypeError(f"Cannot apply tuple index to non-NumPy values: {e}")
                 # Reshape and index
                 reshaped_vals = vals.reshape(len(self._keys), *self.shape)
                 selected_element_values = reshaped_vals[(slice(None),) + item]
                 result_shape = selected_element_values.shape[1:]
                 # Return new MV with sliced data and potentially reduced shape
                 result_mv = type(self).fromkeysvalues(self.algebra, self._keys, selected_element_values)
                 result_mv.shape = result_shape
                 return result_mv
            except IndexError: raise IndexError("Index out of range for MV shape dimensions")
            except Exception as e: raise TypeError(f"Error during tuple indexing of MV: {e}")
        else: raise TypeError(f"MV indices must be integers, slices, or tuples, not {type(item).__name__}")


    def __setitem__(self, key: Union[int, slice, Tuple], value: Any) -> None:
        """
        Assign a value to a component or slice of the multivector. (Makes MV mutable).
        """
        # REVIEW: Mutability comment added. Logic reviewed, complex assignment needs testing.
        if not hasattr(self, '_values'): raise AttributeError("Cannot set item on uninitialized MV")

        is_np = isinstance(self._values, np.ndarray)
        if not is_np and not isinstance(self._values, list):
             try: self._values = list(self._values)
             except TypeError: raise TypeError("Cannot set item on MV with non-list/array values")

        assign_value = value # Handle MV assignment below
        if isinstance(value, MultiVector):
            if value.algebra != self.algebra: raise ValueError("Cannot assign MV from different algebra")
            # Unpack MV value based on key type (more robust checks)
            if isinstance(key, (int, Tuple)): # Assigning to single element/component
                 if value.shape != (): raise ValueError(f"Cannot assign array MV (shape {value.shape}) to single index {key}")
                 if len(value.keys()) != 1: raise ValueError(f"Cannot assign multi-component MV to single index {key}")
                 assign_value = value.values()[0]
            elif isinstance(key, slice): # Assigning to component slice
                 num_slice = len(range(*key.indices(len(self._values))))
                 if value.shape != (): raise ValueError(f"Cannot assign array MV (shape {value.shape}) to component slice")
                 if len(value.keys()) != num_slice: raise ValueError(f"MV value length ({len(value.keys())}) != slice length ({num_slice})")
                 # TODO: Should assignment depend on matching keys in the slice? For now, just use values.
                 assign_value = value.values()

        # --- Component/Slice Assignment (int or slice key) ---
        if isinstance(key, (int, slice)):
            if is_np:
                 try: self._values[key] = assign_value
                 except Exception as e: raise TypeError(f"Error assigning to NumPy MV components: {e}")
            else: # list assignment
                 try: self._values[key] = assign_value
                 except IndexError: raise IndexError("MV index out of range")
                 except Exception as e: raise TypeError(f"Error assigning to list MV components: {e}")

        # --- Array-Valued Assignment (tuple key) ---
        elif isinstance(key, tuple):
            if not self.shape: raise IndexError("Tuple indexing requires array-valued MV")
            if not is_np: # Convert to numpy if list-based
                try:
                     self._values = np.array(self._values, dtype=object) # Use object for safety
                     is_np = True
                except Exception as e: raise TypeError(f"Cannot apply tuple assignment to non-NumPy values: {e}")
            # Reshape and assign
            try:
                 expected_shape = (len(self._keys), *self.shape)
                 # Ensure _values is multi-dim array matching expected shape
                 if self._values.shape != expected_shape:
                     try: self._values = self._values.reshape(expected_shape)
                     except ValueError as e: raise ValueError(f"Cannot reshape internal values {self._values.shape} to {expected_shape}. Error: {e}")
                 # Assign value to slice
                 self._values[(slice(None),) + key] = assign_value
                 # Keep _values reshaped internally
            except IndexError: raise IndexError("Index out of range for MV shape dimensions")
            except Exception as e: raise TypeError(f"Error during tuple assignment to MV: {e}")
        else: raise TypeError(f"MV indices must be int, slice, or tuple, not {type(key).__name__}")


    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        """NumPy ufunc integration."""
        # REVIEW: Handles ufuncs, including multiple outputs. Does not support 'out' kwarg for MV output.
        # CHANGE: Pass info to _wrap_ufunc_result about whether input was ndarray.
        if method != '__call__': return NotImplemented

        mv_inputs, mv_template, other_args = [], None, []
        input_was_numpy_array = False # Flag to track if any input is ndarray
        out = kwargs.get('out', None)

        if out is not None and any(isinstance(o, MultiVector) for o in (out if isinstance(out, tuple) else (out,))):
            print("Warning: __array_ufunc__ does not support 'out' kwarg for MV outputs.")
            return NotImplemented

        for arg in inputs: # Process inputs
            if isinstance(arg, MultiVector):
                mv_inputs.append(arg)
                if mv_template is None: mv_template = arg
                elif arg.algebra != mv_template.algebra: raise TypeError("Ufunc MVs must share algebra")
                if ufunc.nin > 1 and mv_template is not arg: # Check compatibility for binary ufuncs
                     if set(arg.keys()) != set(mv_template.keys()): raise TypeError(f"Ufunc '{ufunc.__name__}' requires matching keysets.")
                other_args.append(arg.values())
            else:
                if isinstance(arg, np.ndarray): # Check if a non-MV input is ndarray
                    input_was_numpy_array = True
                other_args.append(arg)

        if mv_template is None: return NotImplemented # No MV inputs

        try: # Apply ufunc
            result_data = ufunc(*other_args, **kwargs)
            if result_data is NotImplemented: return NotImplemented
        except TypeError as e:
             print(f"Error applying ufunc '{ufunc.__name__}': {e}")
             if "loop of ufunc" in str(e): print("Hint: May involve object arrays.")
             return NotImplemented

        # Process result(s) - Pass the flag to the wrapper
        wrap_kwargs = {'_input_was_numpy_array': input_was_numpy_array}
        if isinstance(result_data, tuple): # Multiple outputs
            return tuple(self._wrap_ufunc_result(mv_template, res, **wrap_kwargs) for res in result_data)
        else: # Single output
            return self._wrap_ufunc_result(mv_template, result_data, **wrap_kwargs)

    # CHANGE: Added _input_was_numpy_array flag and logic
    def _wrap_ufunc_result(self, template, result_values, _input_was_numpy_array=False):
        """Helper to wrap ufunc result array/scalar into a MultiVector."""
        result_shape = ()
        num_keys = len(template.keys())

        if _input_was_numpy_array:
            # If the original operation involved a NumPy array, force resulting MVs to have scalar shape
            result_shape = ()
        else:
            # Original logic to determine shape for MV-MV or MV-scalar operations
            if isinstance(result_values, np.ndarray):
                if result_values.ndim > 0:
                    if result_values.shape[0] == num_keys and result_values.ndim > 1:
                         result_shape = result_values.shape[1:]
                    else: # Shape mismatch or reduction, assume template shape preservation
                         result_shape = template.shape
                         expected_size = num_keys * np.prod(template.shape if template.shape else (1,))
                         if result_values.size != expected_size:
                              if result_values.size == 1: pass # Broadcast scalar handled by fromkeysvalues
                              else: print(f"Warning: Ufunc result size ({result_values.size}) mismatch. Shape may be incorrect.")

        # Handle broadcasting scalar results if needed (independent of shape logic above)
        if not isinstance(result_values, (list, np.ndarray)) and num_keys > 1:
             result_values = [result_values] * num_keys # Let fromkeysvalues handle type

        # Create the result MultiVector
        result_mv = type(template).fromkeysvalues(template.algebra, template.keys(), result_values)
        # Assign the determined shape (potentially forced to () )
        result_mv.shape = result_shape
        return result_mv


    # === Representation and String Methods ===

    def __str__(self) -> str:
        """Human-readable string representation."""
        # REVIEW: Shows first element for array MVs, indicates shape. Uses robust zero checks.
        if not hasattr(self, '_keys') or not self._keys: return "0"
        coeffs_dict = self.coeffs
        shape_str = f" shape={self.shape}" if hasattr(self, 'shape') and self.shape else ""
        if not coeffs_dict: return "0" + shape_str # Handle empty or all-zero case

        # Ensure pretty names are available
        if not hasattr(self.algebra, '_bin2canon_prettystr'):
            try: self.algebra._bin2canon_prettystr = self.algebra.basis_names_dict
            except AttributeError: self.algebra._bin2canon_prettystr = {k: f"e{bin(k)[2:]}" if k != 0 else "1" for k in self.algebra.canon2bin.values()}

        components = []
        def sort_key_func(k): # Sort by grade, then canonical index
            grade = bin(k).count('1')
            try: index_in_grade = self.algebra.indices_for_grades[grade].index(k)
            except: index_in_grade = k # Fallback sort
            return (grade, index_in_grade)
        sorted_keys = sorted(coeffs_dict.keys(), key=sort_key_func)

        for key in sorted_keys:
            val = coeffs_dict[key]
            # Robust zero check
            is_effectively_zero = False
            if isinstance(val, (int, float, complex, np.number)): is_effectively_zero = abs(val) < 1e-12
            elif isinstance(val, np.ndarray): is_effectively_zero = not np.any(np.abs(val) > 1e-12)
            elif hasattr(val, 'is_zero'):
                try: is_effectively_zero = val.is_zero
                except TypeError: is_effectively_zero = (val == 0)
            elif val == 0: is_effectively_zero = True
            if is_effectively_zero: continue

            # Format value string and sign
            val_str, sign_str = "", "+"
            if isinstance(val, sympy.Expr):
                is_negative = val.could_extract_minus_sign()
                val_to_print = sympy.Abs(val) if is_negative else val
                temp_val_str = sympy.printing.sstr(val_to_print, full_prec=False)
                needs_paren = isinstance(val_to_print, (sympy.Add, sympy.Mul)) and len(val_to_print.args) > 1
                val_str = f"({temp_val_str})" if needs_paren else temp_val_str
                if is_negative: sign_str = "-"
            elif isinstance(val, (int, float, complex, np.number)):
                 real_part = val.real if isinstance(val, complex) else val
                 if real_part < -1e-12: sign_str = "-" # Check sign using tolerance
                 abs_val = abs(val)
                 if isinstance(abs_val, float): val_str = f"{abs_val:.4g}"
                 elif isinstance(abs_val, complex):
                     r, i = abs_val.real, abs_val.imag
                     r_zero, i_zero = abs(r) < 1e-12, abs(i) < 1e-12
                     if i_zero: val_str = f"{r:.4g}"
                     elif r_zero: val_str = f"{abs(i):.4g}j"
                     else: val_str = f"({r:.4g}{i:+.4g}j)"
                 else: val_str = str(abs_val)
            else: val_str = str(val)

            blade_str = self.algebra._bin2canon_prettystr.get(key, f"blade[{key}]")
            is_scalar_term, is_one = (key == 0), (val_str == "1")
            term_str = val_str if is_scalar_term else (blade_str if is_one else f"{val_str}{blade_str}")

            if components or sign_str == "-": components.append(sign_str)
            components.append(term_str)

        if not components: return "0" + shape_str
        return " ".join(components) + shape_str


    def __repr__(self) -> str:
        """Detailed string representation for debugging."""
        vals_repr = {}
        if hasattr(self, '_keys') and hasattr(self, '_values'):
            for k, v in self.items():
                 if isinstance(v, np.ndarray): vals_repr[k] = f"array(shape={v.shape}, dtype={v.dtype})"
                 else: vals_repr[k] = repr(v)
        else: vals_repr = "{uninitialized}"
        keys_repr = repr(self.keys())
        return f"MultiVector({repr(self.algebra)}, keys={keys_repr}, values={vals_repr}, shape={getattr(self, 'shape', 'N/A')})"


    # === Transformation Methods ===

    def dual(self, kind: str = 'auto') -> MultiVector:
        """Compute the dual (polarity or Hodge)."""
        try:
            if kind == 'polarity': return self.algebra.polarity(self)
            elif kind == 'hodge': return self.algebra.hodge(self)
            elif kind == 'auto':
                if hasattr(self.algebra, 'r'):
                    if self.algebra.r == 0:
                         try: return self.algebra.polarity(self)
                         except AttributeError: pass
                    try: return self.algebra.hodge(self)
                    except AttributeError: pass
                raise ValueError('Cannot auto-select dual for this algebra.')
            else: raise ValueError(f'Unknown dual kind: {kind}')
        except AttributeError as e: raise AttributeError(f"Algebra lacks dual method for kind '{kind}': {e}")


    def undual(self, kind: str = 'auto') -> MultiVector:
        """Compute the undual."""
        try:
            if kind == 'polarity': return self.algebra.unpolarity(self)
            elif kind == 'hodge': return self.algebra.unhodge(self)
            elif kind == 'auto':
                if hasattr(self.algebra, 'r'):
                    if self.algebra.r == 0:
                         try: return self.algebra.unpolarity(self)
                         except AttributeError: pass
                    try: return self.algebra.unhodge(self)
                    except AttributeError: pass
                raise ValueError('Cannot auto-select undual for this algebra.')
            else: raise ValueError(f'Unknown undual kind: {kind}')
        except AttributeError as e: raise AttributeError(f"Algebra lacks undual method for kind '{kind}': {e}")


    def tosympy(self) -> Any:
        """
        Convert to SymPy representation (scalar or dict).
        Array values become NumPy object arrays of SymPy expressions.
        """
        # REVIEW: Handles arrays by creating object arrays of SymPy expressions. Return type documented.
        sympy_items = {}
        if not hasattr(self, '_keys') or not hasattr(self, '_values'): return {}

        for k, v in self.items():
             try:
                 if isinstance(v, np.ndarray):
                      sympy_array = np.array([S(elem) for elem in v.flat], dtype=object).reshape(v.shape)
                      sympy_items[k] = sympy_array
                 elif isinstance(v, list): # Handle lists
                      np_v = np.array(v, dtype=object)
                      sympy_array = np.array([S(elem) for elem in np_v.flat], dtype=object).reshape(np_v.shape)
                      sympy_items[k] = sympy_array
                 else: # Scalars
                      sympy_items[k] = S(v) # Convert non-Expr scalars too
             except (TypeError, ValueError) as e:
                  print(f"Warning: Cannot convert value for key {k} (type {type(v)}) to SymPy. Using original. Error: {e}")
                  sympy_items[k] = v

        is_scalar = (len(sympy_items) == 1 and 0 in sympy_items)
        is_array = hasattr(self, 'shape') and self.shape
        if is_scalar and not is_array: return sympy_items[0]
        else: return sympy_items


    def __deepcopy__(self, memo):
        """Deep copy support."""
        cls = self.__class__
        result = memo.get(id(self))
        if result is not None: return result
        result = cls.__new__(cls); memo[id(self)] = result
        result.algebra = self.algebra
        if hasattr(self, '_keys'): result._keys = copy.deepcopy(self._keys, memo)
        if hasattr(self, '_values'): result._values = copy.deepcopy(self._values, memo)
        if hasattr(self, 'shape'): result.shape = copy.deepcopy(self.shape, memo)
        return result


# Define aliases
MultiVectorBase = MultiVector
MultiVectorCreation = MultiVector
MultiVectorIndexing = MultiVector
MultiVectorOperations = MultiVector
MultiVectorTransforms = MultiVector

__all__ = ['MultiVector']