"""
Operator Dictionary Module for Geometric Algebra
===============================================
(Docstring unchanged)
"""

from __future__ import annotations

import re
import string
import warnings # Added for warnings
from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, Iterator, List, Optional, Set, Tuple, Type, Union, cast, TYPE_CHECKING # Added TYPE_CHECKING

from sympy import Symbol # Keep Symbol import

# Keep CodegenOutput import here for type hints
from kingdon.codegen import CodegenOutput

# Use TYPE_CHECKING for imports only needed for type hints
if TYPE_CHECKING:
    from kingdon.algebra import Algebra
    from kingdon.multivector import MultiVector

# Define AlgebraError locally if not imported from elsewhere
class AlgebraError(Exception):
    """Custom exception for errors related to Geometric Algebra operations."""
    pass

@dataclass
class OperatorDict(Mapping):
    """
    A dictionary-like object that generates and caches geometric algebra operations.
    (Docstring unchanged)
    """
    name: str
    codegen: Callable[..., Any]
    algebra: 'Algebra'
    operator_dict: Dict[Tuple[Tuple[int, ...], ...], CodegenOutput] = field(
        default_factory=dict, init=False, repr=False)
    codegen_symbolcls: Callable = field(default=Symbol, repr=False) # Retained Symbol

    def __getitem__(self, keys_in: Tuple[Tuple[int, ...], ...]) -> CodegenOutput:
        """ Get or create a compiled operation for the given key pattern. """
        # Input validation: Ensure keys_in is a tuple of tuples of ints
        if not isinstance(keys_in, tuple) or not all(isinstance(k_tuple, tuple) for k_tuple in keys_in):
             raise TypeError(f"OperatorDict key must be a tuple of tuples, got {type(keys_in)}")
        # Further validation could check if ints are valid keys for the algebra if needed

        if keys_in not in self.operator_dict:
            # Moved import inside
            from kingdon.codegen import do_codegen
            from kingdon.multivector import MultiVector # Needed for isinstance check

            # Create symbolic multivectors
            mvs = []
            arg_names = string.ascii_lowercase[:len(keys_in)]
            try:
                for name, keys in zip(arg_names, keys_in):
                    # Use algebra.multivector factory for symbolic creation
                    # Pass name explicitly to trigger symbolic path
                    mvs.append(self.algebra.multivector(name=name, keys=keys))
            except (ValueError, TypeError) as e:
                 # Catch errors during symbolic multivector creation
                 raise AlgebraError(f"Error creating symbolic multivector for operator '{self.name}' with input keys {keys_in}: {e}") from e

            # Generate code for the operation
            try:
                 output = do_codegen(self.codegen, *mvs)
                 self.operator_dict[keys_in] = output
            except NotImplementedError:
                 # Propagate NotImplementedError clearly
                 raise NotImplementedError(f"Symbolic code generation for operator '{self.name}' is not implemented for input keys {keys_in}.")
            except Exception as e:
                 # Catch other codegen errors
                 raise AlgebraError(f"Code generation failed for operator '{self.name}' with input keys {keys_in}: {e}") from e

        return self.operator_dict[keys_in]

    def __iter__(self) -> Iterator[Tuple[Tuple[int, ...], ...]]:
        """Return an iterator over the cached key patterns."""
        return iter(self.operator_dict)

    def __len__(self) -> int:
        """Return the number of cached operations."""
        return len(self.operator_dict)

    def __call__(self, *mvs: 'MultiVector') -> 'MultiVector':
        """ Apply the operation to the given multivectors. """
        from kingdon.multivector import MultiVector # Local import

        # --- Input Validation ---
        if not mvs:
             raise TypeError(f"Operator '{self.name}' requires at least one MultiVector argument.")
        if not all(isinstance(mv, MultiVector) for mv in mvs):
             arg_types = ", ".join(type(arg).__name__ for arg in mvs)
             raise TypeError(f"All arguments to operator '{self.name}' must be MultiVector instances. Got: {arg_types}")
        if not all(mv.algebra == self.algebra for mv in mvs):
             raise ValueError(f"All MultiVector arguments for operator '{self.name}' must belong to the same algebra.")

        # --- Get Compiled Function ---
        keys_in = tuple(mv.keys() for mv in mvs)
        try:
             output = self[keys_in] # Calls __getitem__ which handles codegen/caching
        except (AlgebraError, NotImplementedError, TypeError, ValueError) as e:
             # Propagate errors from __getitem__ with context
             raise type(e)(f"Failed to get/generate function for operator '{self.name}' with input keys {keys_in}: {e}") from e

        # --- Prepare Values and Execute ---
        values_in = tuple(mv.values() for mv in mvs)
        try:
             result_values = output.func(*values_in)
        except Exception as e:
             # Add context to errors during execution of generated code
             raise AlgebraError(f"Runtime error executing generated code for operator '{self.name}' (input keys: {keys_in}, output keys: {output.keys_out}): {e}") from e

        # --- Construct Result ---
        try:
            # Use fromkeysvalues for construction
             return MultiVector.fromkeysvalues(
                 self.algebra, keys=output.keys_out, values=result_values
             )
        except (ValueError, TypeError) as e:
            # Add context to errors during result construction
             raise AlgebraError(f"Error constructing result MultiVector for operator '{self.name}' (output keys: {output.keys_out}): {e}") from e


@dataclass
class UnaryOperatorDict(OperatorDict):
    """ Specialized dictionary for unary operators in Geometric Algebra. """
    # (Docstring unchanged)

    def __getitem__(self, keys_in: Union[Tuple[int, ...], Tuple[Tuple[int, ...]]]) -> CodegenOutput:
        """ Get or create a compiled unary operation for the given key pattern. """
        # Handle single tuple input for convenience during call
        if keys_in and isinstance(keys_in[0], int):
            keys_tuple = (keys_in,) # Ensure it's a tuple of tuples: ((k1, k2,...),)
        elif isinstance(keys_in, tuple) and len(keys_in) == 1 and isinstance(keys_in[0], tuple):
            keys_tuple = keys_in # Already in correct format
        else:
            raise TypeError(f"UnaryOperatorDict key must be a tuple of ints or a tuple containing one tuple of ints, got {type(keys_in)}")

        if keys_tuple not in self.operator_dict:
            # Moved import inside
            from kingdon.codegen import do_codegen
            from kingdon.multivector import MultiVector # For isinstance check

            # Create a symbolic multivector
            try:
                # Use algebra.multivector factory
                 mv = self.algebra.multivector(name='a', keys=keys_tuple[0])
            except (ValueError, TypeError) as e:
                 raise AlgebraError(f"Error creating symbolic multivector for unary operator '{self.name}' with input keys {keys_tuple[0]}: {e}") from e

            # Generate code for the unary operation
            try:
                 output = do_codegen(self.codegen, mv)
                 self.operator_dict[keys_tuple] = output
            except NotImplementedError:
                 raise NotImplementedError(f"Symbolic code generation for unary operator '{self.name}' is not implemented for input keys {keys_tuple[0]}.")
            except Exception as e:
                 raise AlgebraError(f"Code generation failed for unary operator '{self.name}' with input keys {keys_tuple[0]}: {e}") from e

        return self.operator_dict[keys_tuple]

    def __call__(self, mv: 'MultiVector') -> 'MultiVector':
        """ Apply the unary operation to the given multivector. """
        from kingdon.multivector import MultiVector # Local import

        # --- Input Validation ---
        if not isinstance(mv, MultiVector):
             raise TypeError(f"Argument to unary operator '{self.name}' must be a MultiVector instance. Got: {type(mv).__name__}")
        if mv.algebra != self.algebra:
             raise ValueError(f"MultiVector argument for unary operator '{self.name}' does not belong to the correct algebra.")

        # --- Get Compiled Function ---
        keys_in = (mv.keys(),) # Pass keys as tuple of tuples: ((keys,),)
        try:
             output = self[keys_in] # Calls __getitem__
        except (AlgebraError, NotImplementedError, TypeError, ValueError) as e:
             raise type(e)(f"Failed to get/generate function for unary operator '{self.name}' with input keys {keys_in[0]}: {e}") from e

        # --- Execute ---
        try:
             # Call the generated function with the multivector's values
             result_values = output.func(mv.values()) # Unary func takes one arg
        except Exception as e:
             raise AlgebraError(f"Runtime error executing generated code for unary operator '{self.name}' (input keys: {keys_in[0]}, output keys: {output.keys_out}): {e}") from e

        # --- Construct Result ---
        try:
             # Use fromkeysvalues for construction
             return MultiVector.fromkeysvalues(
                 self.algebra, keys=output.keys_out, values=result_values
             )
        except (ValueError, TypeError) as e:
             raise AlgebraError(f"Error constructing result MultiVector for unary operator '{self.name}' (output keys: {output.keys_out}): {e}") from e


@dataclass
class Registry(OperatorDict):
    """ A registry for storing and retrieving compiled custom operators. """
    # (Docstring remains the same)
    # __getitem__ and register_operation could benefit from similar error handling enhancements
    # if custom operators are used extensively and prone to errors.
    # For now, keeping them as they are.

    def __getitem__(self, keys_in: Tuple[Tuple[int, ...], ...]) -> CodegenOutput:
        """ Get or create a compiled custom operation for the given key pattern. """
        # Basic validation
        if not isinstance(keys_in, tuple) or not all(isinstance(k_tuple, tuple) for k_tuple in keys_in):
             raise TypeError(f"Registry key must be a tuple of tuples, got {type(keys_in)}")

        if keys_in not in self.operator_dict:
            # Import moved inside
            from kingdon.codegen import do_codegen

            # Assume codegen func takes MVs like other OperatorDicts
            mvs = []
            arg_names = string.ascii_lowercase[:len(keys_in)]
            try:
                 for name, keys in zip(arg_names, keys_in):
                      mvs.append(self.algebra.multivector(name=name, keys=keys))
            except (ValueError, TypeError) as e:
                 raise AlgebraError(f"Error creating symbolic multivector for custom operator '{self.name}' with input keys {keys_in}: {e}") from e

            try:
                 output = do_codegen(self.codegen, *mvs)
                 self.operator_dict[keys_in] = output
            except NotImplementedError:
                 raise NotImplementedError(f"Symbolic code generation for custom operator '{self.name}' is not implemented for input keys {keys_in}.")
            except Exception as e:
                 raise AlgebraError(f"Code generation failed for custom operator '{self.name}' with input keys {keys_in}: {e}") from e

        return self.operator_dict[keys_in]

    def register_operation(self, operation_name: str,
                         operation_func: Callable[..., Any]) -> None:
        """ Register a custom operation in the Algebra. """
        # (Method unchanged)
        setattr(self.algebra, operation_name, operation_func)


@dataclass
class BladeDict(Mapping):
    """ Dictionary of basis blades for efficient access in Geometric Algebra. """
    # (Docstring remains the same)
    algebra: 'Algebra'
    lazy: bool = False
    blades: Dict[str, 'MultiVector'] = field(
        default_factory=dict, init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        """Initialize the blade dictionary."""
        # (Unchanged)
        if not self.lazy:
            if hasattr(self.algebra, 'canon2bin'):
                for blade_name in self.algebra.canon2bin.keys():
                    try: _ = self[blade_name] # Accessing via __getitem__ populates the cache
                    except (AttributeError, KeyError, ValueError) as e:
                        # Log warning if blade creation fails during eager init
                        warnings.warn(f"Failed to pre-cache blade '{blade_name}': {e}", stacklevel=2)
            else: pass

    def __getitem__(self, basis_blade: str) -> 'MultiVector':
        """Get the multivector representing the specified basis blade."""
        from kingdon.multivector import MultiVector # Local import

        # --- Validation using Algebra._blade2canon ---
        try:
            # This call performs validation and canonicalization
            canonical_name, swaps = self.algebra._blade2canon(basis_blade)
        except ValueError as e:
            # Re-raise as AttributeError or KeyError for consistency with mapping access
            raise KeyError(f"Invalid basis blade name or format: '{basis_blade}'. {e}") from e

        # --- Cache Check & Creation ---
        if canonical_name not in self.blades:
            bin_key = self.algebra.canon2bin.get(canonical_name)
            if bin_key is None:
                # This should not happen if _blade2canon succeeded
                raise KeyError(f"Internal error: Canonical name '{canonical_name}' (from '{basis_blade}') not found in algebra mapping.")

            # Use algebra.multivector factory method
            try:
                # Create the canonical blade (value=1)
                 blade_mv = self.algebra.multivector(keys=(bin_key,), values=[1])
                 self.blades[canonical_name] = blade_mv
            except Exception as e_create:
                 # Handle potential errors during multivector creation
                 raise AlgebraError(f"Failed to create blade '{canonical_name}' via algebra.multivector: {e_create}") from e_create

        # --- Return (with sign adjustment) ---
        blade_mv = self.blades[canonical_name]
        # Return a copy to prevent accidental modification of cached blade? Maybe not needed.
        return blade_mv if swaps % 2 == 0 else -blade_mv # Apply sign flip


    def __getattr__(self, blade: str) -> 'MultiVector':
        """Allow attribute access for basis blades (e.g., alg.blades.e12)."""
        # Try __getitem__ first, which now handles validation & canonicalization
        try:
            return self[blade]
        except KeyError as e:
            # Raise AttributeError if __getitem__ fails (invalid name or other issue)
            raise AttributeError(f"'{type(self).__name__}' object has no attribute '{blade}' or it's not a valid blade name. Original error: {e}") from e

    def __len__(self) -> int:
        """Return the number of cached blades."""
        # Return total number of possible blades, not just cached ones
        return len(self.algebra)

    def __iter__(self) -> Iterator[str]:
        """Return an iterator over all canonical blade names."""
        # Iterate over known canonical names from the algebra
        if hasattr(self.algebra, 'canon2bin'):
             # Sort keys for consistent iteration order
             return iter(sorted(self.algebra.canon2bin.keys(), key=lambda k: self.algebra.canon2bin[k]))
        return iter(sorted(self.blades.keys())) # Fallback to cached blades (less reliable)


    def grade(self, *grades: Union[int, Tuple[int, ...]]) -> Dict[str, 'MultiVector']:
        """Return a dictionary of blades corresponding to the given grade(s)."""
        # (Method unchanged)
        if len(grades) == 1 and isinstance(grades[0], tuple): grades_tuple = grades[0]
        else: grades_tuple = tuple(g for g in grades if isinstance(g, int))
        indices = set()
        for g in grades_tuple: indices.update(self.algebra.indices_for_grades.get((g,), []))
        result = {}
        if hasattr(self.algebra, 'bin2canon'):
            # Sort resulting blades by canonical key order
            sorted_indices = sorted(list(indices), key=self.algebra._canonical_key_sort_func())
            for k in sorted_indices:
                blade_name = self.algebra.bin2canon.get(k)
                if blade_name: result[blade_name] = self[blade_name] # Use __getitem__
        return result

    @property
    def e(self) -> 'MultiVector': # Use quotes
        """Get the scalar unit basis blade (e or e0)."""
        # (Method unchanged)
        try: return self['e']
        except (KeyError, AttributeError):
            try: return self['e0']
            except (KeyError, AttributeError): raise AttributeError("Cannot find scalar blade 'e' or 'e0'.")