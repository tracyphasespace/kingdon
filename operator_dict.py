"""
Operator Dictionary Module for Geometric Algebra
===============================================

This module provides dictionary-like classes that manage operators in Geometric Algebra.
These classes handle code generation, caching, and execution of geometric operations.

Key classes:
- OperatorDict: Base class for binary operators (geometric product, outer product, etc.)
- UnaryOperatorDict: Specialized class for unary operators (reverse, grade involution, etc.)
- Registry: Container for storing custom operators
- BladeDict: Dictionary of basis blades for efficient access

These classes are typically instantiated by the Algebra class and not used directly.
"""

from __future__ import annotations

import re
import string
from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, Iterator, List, Optional, Set, Tuple, Type, Union, cast

from sympy import Symbol

from kingdon.codegen import do_codegen, CodegenOutput


class AlgebraError(Exception):
    """
    Custom exception for errors related to Geometric Algebra operations.
    
    This exception is raised when operations cannot be performed due to
    incompatible types, invalid parameters, or other algebraic constraints.
    """
    pass


@dataclass
class OperatorDict(Mapping):
    """
    A dictionary-like object that generates and caches geometric algebra operations.
    
    This class handles binary operators like geometric product, outer product, etc.
    It performs code generation for operations between multivectors with specific
    key patterns, and caches the results for efficiency.
    
    Attributes:
        name: Name of the operator (e.g., "gp", "op", "ip").
        codegen: Function that generates code for the operation.
        algebra: The Geometric Algebra instance this operator belongs to.
        operator_dict: Cache for storing generated operations.
        codegen_symbolcls: Symbol class used during code generation.
    """
    name: str
    codegen: Callable[..., Any]
    algebra: Any  # Type should be 'Algebra' but avoiding circular import
    operator_dict: Dict[Tuple[Tuple[int, ...], ...], CodegenOutput] = field(
        default_factory=dict, init=False, repr=False)
    codegen_symbolcls: Callable = field(default=Symbol, repr=False)

    def __getitem__(self, keys_in: Tuple[Tuple[int, ...], ...]) -> CodegenOutput:
        """
        Get or create a compiled operation for the given key pattern.
        
        Args:
            keys_in: Tuple of key tuples, one for each input multivector.
                    For example, ((1, 2), (4, 8)) for two multivectors with 
                    those specific keys.
        
        Returns:
            A CodegenOutput containing the output keys and callable function.
        """
        if keys_in not in self.operator_dict:
            # Create symbolic multivectors for each set of keys
            mvs = [
                self.algebra.multivector(
                    name=name, keys=keys, symbolcls=self.codegen_symbolcls
                )
                for name, keys in zip(string.ascii_lowercase, keys_in)
            ]
            # Generate code for the operation
            output = do_codegen(self.codegen, *mvs)
            self.operator_dict[keys_in] = output
        return self.operator_dict[keys_in]

    def __iter__(self) -> Iterator[Tuple[Tuple[int, ...], ...]]:
        """Return an iterator over the cached key patterns."""
        return iter(self.operator_dict)

    def __len__(self) -> int:
        """Return the number of cached operations."""
        return len(self.operator_dict)

    def __call__(self, *mvs: Any) -> Any:
        """
        Apply the operation to the given multivectors.
        
        This method retrieves or generates the appropriate compiled function
        for the multivectors' key patterns, then applies it to their values.
        
        Args:
            *mvs: Multivector objects to operate on.
        
        Returns:
            A new multivector containing the result of the operation.
            
        Example:
            >>> result = gp_dict(vector1, vector2)  # Compute geometric product
        """
        from kingdon.multivector import MultiVector  # Local import to avoid circular imports
        
        # Get the key pattern for each multivector
        keys_in = tuple(mv.keys() for mv in mvs)
        
        # Get or generate the compiled operation
        output = self[keys_in]
        
        # Call the generated function with each multivector's values
        values_in = tuple(mv.values() for mv in mvs)
        result_values = output.func(*values_in)
        
        # Create a new multivector with the result
        return MultiVector.fromkeysvalues(
            self.algebra, keys=output.keys_out, values=result_values
        )


@dataclass
class UnaryOperatorDict(OperatorDict):
    """
    Specialized dictionary for unary operators in Geometric Algebra.
    
    This class handles operations that act on a single multivector, such as
    reverse, grade involution, Clifford conjugate, square root, etc.
    
    Attributes:
        Same as OperatorDict, but operations take only one multivector.
    """
    def __getitem__(self, keys_in: Union[Tuple[int, ...], Tuple[Tuple[int, ...]]]) -> CodegenOutput:
        """
        Get or create a compiled unary operation for the given key pattern.
        
        Args:
            keys_in: Either a tuple of keys, or a tuple containing a tuple of keys.
                    For example, (1, 2, 4) or ((1, 2, 4),).
        
        Returns:
            A CodegenOutput containing the output keys and callable function.
        """
        # Handle different input formats for convenience
        if keys_in and isinstance(keys_in[0], int):
            keys_in = (keys_in,)  # type: ignore
            
        # Store as a tuple containing a tuple for consistency with mapping
        keys_tuple = cast(Tuple[Tuple[int, ...]], keys_in)
        
        if keys_tuple not in self.operator_dict:
            # Create a symbolic multivector for the key pattern
            mv = self.algebra.multivector(
                name='a', keys=keys_tuple[0], symbolcls=self.codegen_symbolcls
            )
            # Generate code for the unary operation
            output = do_codegen(self.codegen, mv)
            self.operator_dict[keys_tuple] = output
            
        return self.operator_dict[keys_tuple]

    def __call__(self, mv: Any) -> Any:
        """
        Apply the unary operation to the given multivector.
        
        Args:
            mv: The multivector to operate on.
        
        Returns:
            A new multivector containing the result of the operation.
            
        Example:
            >>> result = reverse_dict(bivector)  # Compute reverse of bivector
        """
        from kingdon.multivector import MultiVector  # Local import to avoid circular imports
        
        # Get or generate the compiled operation
        output = self[(mv.keys(),)]
        
        # Call the generated function with the multivector's values
        result_values = output.func(mv.values())
        
        # Create a new multivector with the result
        return MultiVector.fromkeysvalues(
            self.algebra, keys=output.keys_out, values=result_values
        )


@dataclass
class Registry(OperatorDict):
    """
    A registry for storing and retrieving compiled custom operators.
    
    This class allows extending the Algebra with user-defined operations.
    
    Attributes:
        Same as OperatorDict, plus additional methods for registering
        custom operations.
    """
    def __getitem__(self, keys_in: Tuple[Tuple[int, ...], ...]) -> CodegenOutput:
        """
        Get or create a compiled custom operation for the given key pattern.
        
        Args:
            keys_in: Tuple of key tuples, one for each input multivector.
        
        Returns:
            A CodegenOutput containing the output keys and callable function.
        """
        if keys_in not in self.operator_dict:
            # For Registry, we use do_codegen with keys_in directly
            output = do_codegen(self.codegen, *keys_in)
            self.operator_dict[keys_in] = output
        return self.operator_dict[keys_in]
    
    def register_operation(self, operation_name: str, 
                         operation_func: Callable[..., Any]) -> None:
        """
        Register a custom operation in the Algebra.
        
        Args:
            operation_name: Name for the custom operation.
            operation_func: Function implementing the operation.
            
        Example:
            >>> # Define a custom operation
            >>> def my_custom_op(x, y):
            ...     return x * y + x | y
            >>>
            >>> # Register it in the algebra
            >>> registry.register_operation("custom_op", my_custom_op)
            >>> 
            >>> # Now it can be used like other operators
            >>> result = algebra.custom_op(vector1, vector2)
        """
        # Implementation would integrate the operation with the algebra
        # This is a placeholder for the actual implementation
        setattr(self.algebra, operation_name, operation_func)


@dataclass
class BladeDict(Mapping):
    """
    Dictionary of basis blades for efficient access in Geometric Algebra.
    
    This class lazily creates and caches multivectors representing basis blades
    when they are first accessed. This improves performance, especially in
    higher-dimensional algebras.
    
    Attributes:
        algebra: The Geometric Algebra instance.
        lazy: If True, blades are created only when accessed. Otherwise,
              all blades are pre-computed during initialization.
        blades: Cache storing created blade multivectors.
    
    Usage:
        >>> blade_dict = BladeDict(algebra=alg, lazy=True)
        >>> blade = blade_dict['e12']   # Retrieves the blade for e12
        >>> blade2 = blade_dict.e12     # Also accessible as an attribute
    """
    algebra: Any  # Type should be 'Algebra' but avoiding circular import
    lazy: bool = False
    blades: Dict[str, Any] = field(
        default_factory=dict, init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        """
        Initialize the blade dictionary.
        
        If not lazy, pre-compute all basis blades.
        """
        if not self.lazy:
            # Pre-load all blades if not lazy
            for blade in self.algebra.canon2bin:
                _ = self[blade]

    def __getitem__(self, basis_blade: str) -> Any:
        """
        Get the multivector representing the specified basis blade.
        
        Args:
            basis_blade: String representing the basis blade (e.g., 'e12', 'e0', 'e123').
        
        Returns:
            Multivector representing the basis blade.
            
        Raises:
            AttributeError: If the basis blade string is invalid.
        """
        from kingdon.multivector import MultiVector  # Local import to avoid circular imports
        
        # Validate the basis blade string format
        if not re.match(r'^e[0-9a-fA-F]*$', basis_blade):
            raise AttributeError(f"'{basis_blade}' is not a valid basis blade.")
            
        # Special handling for scalar basis (e0 or e)
        if basis_blade == "e0" or basis_blade == "e":
            canonical, swaps = basis_blade, 0
        else:
            # Check if _blade2canon is available
            if hasattr(self.algebra, '_blade2canon'):
                # Convert to canonical form and determine if a sign swap is needed
                canonical, swaps = self.algebra._blade2canon(basis_blade)
            else:
                # Simple fallback: use the blade as-is, no swaps
                canonical, swaps = basis_blade, 0
        
        # Create the blade if not already cached
        if canonical not in self.blades:
            # Get blade's binary key
            if canonical in self.algebra.canon2bin:
                bin_blade = self.algebra.canon2bin[canonical]
            else:
                # Extract indices and convert to binary key
                indices = [int(c) for c in canonical[1:]] if len(canonical) > 1 else [0]
                bin_blade = sum(1 << (i - self.algebra.start_index) for i in indices)
            
            # Create the multivector for this blade
            if hasattr(self.algebra, 'multivector'):
                # Use the algebra's multivector method if available
                self.blades[canonical] = self.algebra.multivector(
                    values=[1], keys=(bin_blade,)
                )
            else:
                # Direct construction as fallback
                mv = object.__new__(MultiVector)
                mv.algebra = self.algebra
                mv._keys = (bin_blade,)
                mv._values = [1]
                mv.shape = ()
                self.blades[canonical] = mv
                
        blade_mv = self.blades[canonical]
        
        # Apply sign flip if needed due to reordering
        return blade_mv if swaps % 2 == 0 else -blade_mv

    def __getattr__(self, blade: str) -> Any:
        """
        Allow attribute access for basis blades.
        
        This enables accessing blades as attributes, e.g., blades.e12
        instead of blades['e12'].
        
        Args:
            blade: Name of the basis blade.
            
        Returns:
            Multivector representing the basis blade.
            
        Raises:
            AttributeError: If the attribute is not a valid basis blade.
        """
        if blade.startswith('e'):
            return self[blade]
        raise AttributeError(f"'{self.__class__.__name__}' object has no attribute '{blade}'")

    def __len__(self) -> int:
        """Return the number of cached blades."""
        return len(self.blades)

    def __iter__(self) -> Iterator[str]:
        """Return an iterator over the cached blade names."""
        return iter(self.blades)

    def grade(self, *grades: Union[int, Tuple[int, ...]]) -> Dict[str, Any]:
        """
        Return a dictionary of blades corresponding to the given grade(s).
        
        Args:
            *grades: The grade(s) to filter on. Can be individual integers
                    or a tuple of integers.
        
        Returns:
            Dictionary mapping canonical blade names to their multivectors.
            
        Example:
            >>> # Get all bivectors (grade 2)
            >>> bivectors = blades.grade(2)
            >>> 
            >>> # Get all scalars and vectors (grades 0 and 1)
            >>> scalars_and_vectors = blades.grade(0, 1)
        """
        from kingdon.multivector import MultiVector  # Local import to avoid circular imports
        
        # Handle the case where grades is a tuple within a tuple
        if len(grades) == 1 and isinstance(grades[0], tuple):
            grades = grades[0]
            
        # Convert to tuple of integers
        grade_tuple = tuple(
            g for g in grades if isinstance(g, int)
        )
        
        # Get the indices for the specified grades
        indices = self.algebra.indices_for_grades.get(grade_tuple, [])
        
        # Build and return the dictionary of blades
        result = {}
        for k in indices:
            if k in self.algebra.bin2canon:
                blade = self.algebra.bin2canon[k]
                result[blade] = self[blade]
                
        return result
    
    @property
    def e(self) -> Any:
        """
        Get the scalar unit basis blade (e or e0).
        
        Returns:
            Multivector representing the scalar unit.
        """
        # Try 'e' first, then 'e0' as fallback
        try:
            return self['e']
        except (KeyError, AttributeError):
            return self['e0']