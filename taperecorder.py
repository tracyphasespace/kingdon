"""
TapeRecorder Module for Geometric Algebra
=========================================

This module provides the TapeRecorder class, which creates a symbolic representation
of geometric algebra operations for delayed evaluation and code generation.

Unlike MultiVector, which computes values immediately, TapeRecorder creates a
string representation of operations that can be executed later or converted to code.
This enables efficient computation pipelines and code generation for performance-critical
applications.

Example:
    >>> from kingdon.algebra import Algebra
    >>> from kingdon.taperecorder import TapeRecorder
    >>> 
    >>> alg = Algebra(p=3, q=0, r=0)  # 3D Euclidean GA
    >>> # Create a recorder for a vector expression
    >>> x = TapeRecorder(alg, "x", (1, 2, 4))  # Keys for e1, e2, e3
    >>> y = TapeRecorder(alg, "y", (1, 2, 4))
    >>> 
    >>> # Record operations without evaluating them
    >>> z = x * y  # Creates expression representing geometric product
    >>> print(z.expr)
    # Output shows function call: gp_1_2_4_x_1_2_4(x, y)
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from functools import cached_property, partial, partialmethod
from typing import Any, Callable, Dict, List, Optional, Sequence, Set, Tuple, Union, TYPE_CHECKING

# Use TYPE_CHECKING to avoid circular imports
if TYPE_CHECKING:
    from kingdon.algebra import Algebra


@dataclass(init=False)
class TapeRecorder:
    """
    Records geometric algebra operations as symbolic expressions for delayed evaluation.
    
    TapeRecorder allows for creating computation pipelines without immediate execution,
    enabling code generation and optimization. It mimics the MultiVector interface
    but records operations rather than computing them.
    
    Attributes:
        algebra: The Geometric Algebra instance this recorder operates in
        expr: String expression representing the computation
        _keys: Tuple of keys representing the basis blades in the expression
        
    Example:
        >>> x = TapeRecorder(algebra, "x", (1, 2, 4))  # e1, e2, e3 vector
        >>> y = TapeRecorder(algebra, "y", (1, 2, 4))
        >>> z = x * y  # Records the geometric product
        >>> print(z.expr)
        gp_1_2_4_x_1_2_4(x, y)
    """
    algebra: Algebra
    expr: str
    _keys: Tuple[int, ...] = field(default_factory=tuple)

    def __new__(cls, algebra: Any, expr: str, keys: Tuple[int, ...]) -> TapeRecorder:
        """
        Create a new TapeRecorder instance.
        
        Args:
            algebra: The Geometric Algebra instance
            expr: String expression representing the computation
            keys: Tuple of keys representing the basis blades
            
        Returns:
            A new TapeRecorder instance
        """
        obj = object.__new__(cls)
        obj.algebra = algebra
        obj.expr = expr
        obj._keys = keys
        return obj

    def keys(self) -> Tuple[int, ...]:
        """
        Get the basis blade keys in this expression.
        
        Returns:
            Tuple of integer keys representing the basis blades
        """
        return self._keys

    @cached_property
    def type_number(self) -> int:
        """
        Compute a unique type number based on which basis blades are present.
        
        Returns:
            Integer representing the type of this expression
        """
        # Create a binary string with 1s for present keys and 0s for absent keys
        return int(''.join('1' if i in self.keys() else '0' 
                           for i in reversed(self.algebra.canon2bin.values())), 2)

    def __getattr__(self, basis_blade: str) -> TapeRecorder:
        """
        Extract a component corresponding to a specific basis blade.
        
        Args:
            basis_blade: Name of the basis blade (e.g., 'e12')
            
        Returns:
            A new TapeRecorder representing just that component
            
        Raises:
            AttributeError: If the string doesn't match a basis blade pattern
        """
        if not re.match(r'^e[0-9a-fA-F]*$', basis_blade):
            raise AttributeError(
                f'{self.__class__.__name__} object has no attribute or basis blade {basis_blade}')
            
        # Check if basis_blade is in the algebra
        if basis_blade not in self.algebra.canon2bin:
            return self.__class__(
                algebra=self.algebra,
                expr="(0,)",
                keys=(0,)
            )
            
        # Try to find the index of the basis blade in our keys
        try:
            idx = self.keys().index(self.algebra.canon2bin[basis_blade])
        except ValueError:
            # Basis blade not present in this expression
            return self.__class__(
                algebra=self.algebra,
                expr="(0,)",
                keys=(0,)
            )
        else:
            # Extract the component corresponding to the basis blade
            return self.__class__(
                algebra=self.algebra,
                expr=f"({self.expr}[{idx}],)",
                keys=(self.keys()[idx],)
            )

    def grade(self, *grades: Union[int, Tuple[int, ...]]) -> TapeRecorder:
        """
        Extract components of specific grades.
        
        Args:
            *grades: One or more grades to extract, or a tuple of grades
            
        Returns:
            A new TapeRecorder containing only the specified grades
            
        Example:
            >>> x = TapeRecorder(algebra, "x", (0, 1, 2, 4, 8))
            >>> # Extract scalar and vector parts (grades 0 and 1)
            >>> x_sv = x.grade(0, 1)
        """
        # Handle case where grades is passed as a tuple
        if len(grades) == 1 and isinstance(grades[0], tuple):
            grades = grades[0]

        # Get the indices of basis blades with the specified grades
        basis_blades = self.algebra.indices_for_grades[grades]
        indices_keys = [(idx, k) for idx, k in enumerate(self.keys()) if k in basis_blades]
        
        # Extract the matching indices and keys
        indices, keys = zip(*indices_keys) if indices_keys else (tuple(), tuple())
        expr = f"[{self.expr}[idx] for idx in {indices}]"
        
        return self.__class__(
            algebra=self.algebra,
            expr=expr,
            keys=keys,
        )

    def __str__(self) -> str:
        """
        Return the string representation of this expression.
        
        Returns:
            The expression string
        """
        return self.expr

    def binary_operator(self, other: Any, operator: str) -> TapeRecorder:
        """
        Apply a binary operator to this expression and another.
        
        Args:
            other: Another TapeRecorder or a scalar
            operator: Name of the operator to apply
            
        Returns:
            A new TapeRecorder representing the operation
        """
        if not isinstance(other, self.__class__):
            # Assume other is a scalar
            keys_out, func = getattr(self.algebra, operator)[self.keys(), (0,)]
            expr = f'{func.__name__}({self.expr}, ({other},))'
        else:
            # Both are TapeRecorders
            keys_out, func = getattr(self.algebra, operator)[self.keys(), other.keys()]
            expr = f'{func.__name__}({self.expr}, {other.expr})'
            
        return self.__class__(algebra=self.algebra, expr=expr, keys=keys_out)

    def unary_operator(self, operator: str) -> TapeRecorder:
        """
        Apply a unary operator to this expression.
        
        Args:
            operator: Name of the operator to apply
            
        Returns:
            A new TapeRecorder representing the operation
        """
        keys_out, func = getattr(self.algebra, operator)[self.keys()]
        expr = f'{func.__name__}({self.expr})'
        return self.__class__(algebra=self.algebra, expr=expr, keys=keys_out)

    # Binary operators
    gp = __mul__ = __rmul__ = partialmethod(binary_operator, operator='gp')
    sw = __rshift__ = partialmethod(binary_operator, operator='sw')
    cp = partialmethod(binary_operator, operator='cp')
    acp = partialmethod(binary_operator, operator='acp')
    ip = __or__ = partialmethod(binary_operator, operator='ip')
    sp = partialmethod(binary_operator, operator='sp')
    lc = partialmethod(binary_operator, operator='lc')
    rc = partialmethod(binary_operator, operator='rc')
    op = __xor__ = __rxor__ = partialmethod(binary_operator, operator='op')
    rp = __and__ = partialmethod(binary_operator, operator='rp')
    proj = __matmul__ = partialmethod(binary_operator, operator='proj')
    add = __add__ = __radd__ = partialmethod(binary_operator, operator='add')
    sub = __sub__ = partialmethod(binary_operator, operator='sub')
    
    def __rsub__(self, other: Any) -> TapeRecorder:
        """
        Subtract this expression from something else.
        
        Args:
            other: Value to subtract from (typically a scalar)
            
        Returns:
            A new TapeRecorder representing the operation
        """
        # Convert other to a TapeRecorder if it's not already
        if not isinstance(other, TapeRecorder):
            # Create a scalar TapeRecorder
            other_recorder = TapeRecorder(
                self.algebra, 
                f"({other},)", 
                (0,)  # Scalar has key 0
            )
        else:
            other_recorder = other
            
        return other_recorder + (-self)
    
    __truediv__ = div = partialmethod(binary_operator, operator='div')

    def __pow__(self, power: int, modulo: Any = None) -> TapeRecorder:
        """
        Raise this expression to a power.
        
        Args:
            power: Integer exponent
            modulo: Not used, included for compatibility
            
        Returns:
            A new TapeRecorder representing the operation
            
        Raises:
            ValueError: If power is negative
        """
        if power < 0:
            raise ValueError("Negative powers not supported directly. Use inv() first.")
            
        if power == 0:
            # Return a scalar 1
            return self.__class__(
                self.algebra, 
                expr='(1,)', 
                keys=(0,)
            )

        # For powers > 0, multiply the expression by itself repeatedly
        res = self
        for i in range(1, power):
            res = res.gp(self)
            
        return res

    # Unary operators
    inv = partialmethod(unary_operator, operator='inv')
    neg = __neg__ = partialmethod(unary_operator, operator='neg')
    reverse = __invert__ = partialmethod(unary_operator, operator='reverse')
    involute = partialmethod(unary_operator, operator='involute')
    conjugate = partialmethod(unary_operator, operator='conjugate')
    sqrt = partialmethod(unary_operator, operator='sqrt')
    polarity = partialmethod(unary_operator, operator='polarity')
    unpolarity = partialmethod(unary_operator, operator='unpolarity')
    hodge = partialmethod(unary_operator, operator='hodge')
    unhodge = partialmethod(unary_operator, operator='unhodge')
    normsq = partialmethod(unary_operator, operator='normsq')
    outerexp = partialmethod(unary_operator, operator='outerexp')
    outersin = partialmethod(unary_operator, operator='outersin')
    outercos = partialmethod(unary_operator, operator='outercos')
    outertan = partialmethod(unary_operator, operator='outertan')

    def dual(self, kind: str = 'auto') -> TapeRecorder:
        """
        Compute the dual of this expression.
        
        Args:
            kind: Type of dual to compute ('polarity', 'hodge', or 'auto')
            
        Returns:
            A new TapeRecorder representing the dual
            
        Raises:
            ValueError: If the requested dual type is not available
            Exception: If 'auto' can't determine a suitable dual
        """
        if kind == 'polarity' or (kind == 'auto' and self.algebra.r == 0):
            return self.polarity()
        elif kind == 'hodge' or (kind == 'auto' and self.algebra.r == 1):
            return self.hodge()
        elif kind == 'auto':
            raise Exception('Cannot select a suitable dual in auto mode for this algebra.')
        else:
            raise ValueError(f'No dual found for kind={kind}.')

    def undual(self, kind: str = 'auto') -> TapeRecorder:
        """
        Compute the undual (inverse of dual) of this expression.
        
        Args:
            kind: Type of undual to compute ('polarity', 'hodge', or 'auto')
            
        Returns:
            A new TapeRecorder representing the undual
            
        Raises:
            ValueError: If the requested undual type is not available
            Exception: If 'auto' can't determine a suitable undual
        """
        if kind == 'polarity' or (kind == 'auto' and self.algebra.r == 0):
            return self.unpolarity()
        elif kind == 'hodge' or (kind == 'auto' and self.algebra.r == 1):
            return self.unhodge()
        elif kind == 'auto':
            raise Exception('Cannot select a suitable undual in auto mode for this algebra.')
        else:
            raise ValueError(f'No undual found for kind={kind}.')

    def norm(self) -> TapeRecorder:
        """
        Compute the norm (square root of the squared norm) of this expression.
        
        Returns:
            A new TapeRecorder representing the norm
        """
        normsq = self.normsq()
        return normsq.sqrt()

    def normalized(self) -> TapeRecorder:
        """
        Compute the normalized version of this expression.
        
        Returns:
            A new TapeRecorder representing the normalized expression
        """
        return self / self.norm()
    
    def compile(self) -> Callable:
        """
        Compile this expression into an executable function.
        
        This function creates Python code from the tape recording and compiles
        it into a callable function.
        
        Returns:
            A callable function that executes the recorded operations
        """
        # This is a placeholder for actual compilation logic
        # In a real implementation, this would generate code, compile it,
        # and return a callable function
        
        # Example implementation:
        # import types
        # code = f"def _compiled_func(*args):\n    return {self.expr}"
        # namespace = {}
        # exec(code, namespace)
        # return namespace['_compiled_func']
        
        raise NotImplementedError("Compilation not yet implemented")
    
    @property
    def e(self) -> int:
        """
        Get the scalar (e0) component of this expression.
        
        Returns:
            The scalar component as a string expression
            
        Raises:
            AttributeError: If the scalar component doesn't exist
        """
        return getattr(self, "e0")