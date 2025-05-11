"""
TapeRecorder Module for Geometric Algebra
=========================================
(Docstring unchanged)
"""

from __future__ import annotations

import re
import warnings # Added for warnings
from dataclasses import dataclass, field
from functools import cached_property, partial, partialmethod
from typing import Any, Callable, Dict, List, Optional, Sequence, Set, Tuple, Union, TYPE_CHECKING

# Use TYPE_CHECKING to avoid circular imports
if TYPE_CHECKING:
    from kingdon.algebra import Algebra, AlgebraError
    from kingdon.operator_dict import OperatorDict, UnaryOperatorDict # For type hints

# Ensure AlgebraError is defined or imported
try:
    from kingdon.operator_dict import AlgebraError
except ImportError:
    class AlgebraError(Exception):
        """Custom exception for errors related to Geometric Algebra operations."""
        pass


@dataclass(init=False)
class TapeRecorder:
    """
    Records geometric algebra operations as symbolic expressions for delayed evaluation.
    (Docstring unchanged)

    Attributes:
        algebra: The Geometric Algebra instance this recorder operates in.
        expr: String expression representing the computation graph node.
        _keys: Tuple of integer keys representing the basis blades potentially
               present in the result of this expression node.
    """
    algebra: 'Algebra'
    expr: str
    _keys: Tuple[int, ...] = field(default_factory=tuple)

    # Using __init__ instead of __new__ for standard dataclass initialization
    def __init__(self, algebra: 'Algebra', expr: str, keys: Tuple[int, ...]):
        """
        Initialize a TapeRecorder instance.

        Args:
            algebra: The Geometric Algebra instance.
            expr: String expression for this node.
            keys: Tuple of potential output keys.

        Raises:
            TypeError: If inputs have incorrect types.
            ValueError: If keys are invalid for the algebra.
        """
        if not hasattr(algebra, 'd'): # Basic check for Algebra-like object
             raise TypeError("`algebra` argument must be a valid kingdon Algebra instance.")
        if not isinstance(expr, str):
             raise TypeError(f"`expr` must be a string, got {type(expr).__name__}")
        if not isinstance(keys, tuple) or not all(isinstance(k, int) for k in keys):
             raise TypeError(f"`keys` must be a tuple of integers, got {type(keys).__name__}")

        # Validate keys against algebra dimension
        max_key = (1 << algebra.d) - 1
        if not all(0 <= k <= max_key for k in keys):
             invalid_keys = [k for k in keys if not (0 <= k <= max_key)]
             raise ValueError(f"Keys {invalid_keys} are out of range [0, {max_key}] for algebra dimension {algebra.d}.")

        self.algebra = algebra
        self.expr = expr
        # Sort keys for internal consistency? Matches MultiVector.
        self._keys = tuple(sorted(list(set(keys)), key=algebra._canonical_key_sort_func()))


    def keys(self) -> Tuple[int, ...]:
        """ Get the basis blade keys potentially present in this expression. """
        return self._keys

    @cached_property
    def type_number(self) -> int:
        """ Compute a unique type number based on potential basis blades. """
        # (Method unchanged)
        present_keys = set(self.keys()) # Use the stored keys
        # Ensure canon2bin is initialized if needed
        if not hasattr(self.algebra, 'canon2bin') or not self.algebra.canon2bin:
            self.algebra._setup_basis_mappings()

        canon_keys_ordered = list(self.algebra.canon2bin.values())
        # Handle algebras where canon2bin might be smaller than 2**d?
        # Safer to iterate 0 to 2**d - 1 if canon2bin might be incomplete
        # max_key = (1 << self.algebra.d)
        # bitstr = "".join("1" if key in present_keys else "0" for key in range(max_key - 1, -1, -1))

        # Original logic relying on canon2bin values order:
        bitstr = "".join("1" if key in present_keys else "0"
                           for key in reversed(canon_keys_ordered))
        return int(bitstr, 2) if bitstr else 0


    def __getattr__(self, basis_blade_name: str) -> TapeRecorder:
        """
        Extract a component corresponding to a specific basis blade name.
        Returns a *scalar* TapeRecorder representing that component's expression,
        or a zero scalar TapeRecorder if the blade is not present.
        """
        # Validate basis_blade_name using algebra's method
        try:
            canonical_name, _ = self.algebra._blade2canon(basis_blade_name)
            target_key = self.algebra.canon2bin.get(canonical_name)
            if target_key is None:
                # Should not happen if _blade2canon succeeded
                raise KeyError(f"Internal error: Cannot find key for canonical name '{canonical_name}'.")
        except (ValueError, KeyError) as e:
             # Raise AttributeError for consistency with object attribute access
             raise AttributeError(f"'{type(self).__name__}' object has no attribute '{basis_blade_name}' or it's not a valid blade name. Error: {e}") from e

        # Check if the target key is among the keys potentially produced by this expression
        try:
            # Find the index corresponding to the target_key within self._keys
            idx = self.keys().index(target_key)
            # Create expression to extract this specific component
            # Assumes self.expr evaluates to a sequence (list/tuple)
            component_expr = f"{self.expr}[{idx}]"
            # The result represents a single scalar component
            return TapeRecorder(algebra=self.algebra, expr=f"({component_expr},)", keys=(0,))
        except ValueError:
            # Basis blade is valid but not present in this expression's keys
            # Return a zero scalar TapeRecorder
            return TapeRecorder(algebra=self.algebra, expr="(0,)", keys=(0,))


    def grade(self, *grades: Union[int, Tuple[int, ...]]) -> TapeRecorder:
        """ Extract components of specific grades. """
        # Handle single tuple argument e.g. grade((1,3))
        if len(grades) == 1 and isinstance(grades[0], (tuple, list, set)):
            grades_set = set(grades[0])
        else:
            grades_set = set()
            for g in grades:
                 if not isinstance(g, int): raise TypeError(f"Grades must be integers, got {type(g).__name__}")
                 grades_set.add(g)

        # Validate grades
        if not all(0 <= g <= self.algebra.d for g in grades_set):
             invalid_grades = [g for g in grades_set if not (0 <= g <= self.algebra.d)]
             raise ValueError(f"Invalid grades {invalid_grades}. Grades must be between 0 and {self.algebra.d}.")

        # Find which of *self's* keys match the requested grades
        matching_indices_keys = []
        current_keys = self.keys() # Keys potentially output by self.expr
        for idx, k in enumerate(current_keys):
             if bin(k).count('1') in grades_set:
                 matching_indices_keys.append((idx, k))

        if not matching_indices_keys:
            # Return zero scalar if no components match
            return TapeRecorder(algebra=self.algebra, expr="(0,)", keys=(0,))

        # Extract the matching indices and keys
        indices, keys_out = zip(*matching_indices_keys)
        # Sort output keys canonically
        keys_out_sorted = tuple(sorted(keys_out, key=self.algebra._canonical_key_sort_func()))

        # Build expression to select and potentially reorder components
        # If keys_out == keys_out_sorted, we just need to select
        if keys_out == keys_out_sorted:
             # Simple selection using list comprehension of indices
             if len(indices) == 1: # Single component selected
                  expr = f"({self.expr}[{indices[0]}],)" # Wrap in tuple expr
             else:
                  expr = f"[{self.expr}[idx] for idx in {indices}]"
        else:
             # Need to reorder. Build a dict mapping original index to value, then extract sorted.
             # This is complex to represent purely as a string expression without intermediate vars.
             # For simplicity, we might just select and accept the non-canonical order from self.expr,
             # or implement a helper function in numspace for reordering.
             # Let's keep the selection simple for now, the keys_out_sorted are the important part.
             warnings.warn("TapeRecorder.grade() output expression may not reflect canonical key order.", stacklevel=2)
             if len(indices) == 1:
                  expr = f"({self.expr}[{indices[0]}],)"
             else:
                  expr = f"[{self.expr}[idx] for idx in {indices}]"
             # Use the unsorted keys corresponding to the selected expression order
             keys_out_final = keys_out
             # Or raise NotImplementedError("Reordering needed for TapeRecorder.grade, not implemented.")

        return TapeRecorder(
            algebra=self.algebra,
            expr=expr,
            keys=keys_out_sorted, # Report the canonical keys for the result type
        )

    def __str__(self) -> str:
        """ Return the string representation of this expression node. """
        return self.expr

    # --- Operator Implementation ---

    def _get_operator_output(self, operator_name: str, *other_keys_tuples: Tuple[Tuple[int,...]]):
        """ Safely get operator CodegenOutput, raising informative errors. """
        try:
            op_dict = getattr(self.algebra, operator_name)
            # Construct the key for the operator dictionary lookup
            full_keys_tuple = (self.keys(),) + other_keys_tuples
            return op_dict[full_keys_tuple]
        except AttributeError:
            raise NotImplementedError(f"Operator '{operator_name}' not defined in algebra {self.algebra}.")
        except (TypeError, ValueError, KeyError, NotImplementedError, AlgebraError) as e:
            # Catch errors from OperatorDict.__getitem__ (e.g., codegen failure)
            input_keys_repr = (self.keys(),) + other_keys_tuples
            raise AlgebraError(f"Failed to get generated function for '{operator_name}' with keys {input_keys_repr}: {e}") from e


    def binary_operator(self, other: Any, operator: str) -> TapeRecorder:
        """ Apply a binary operator: self OP other """
        other_expr_str: str
        other_keys_tuple: Tuple[int, ...]

        if isinstance(other, TapeRecorder):
            if self.algebra != other.algebra:
                 raise ValueError(f"Cannot apply operator '{operator}' between TapeRecorders from different algebras.")
            other_expr_str = other.expr
            other_keys_tuple = other.keys()
        elif isinstance(other, (int, float, complex)): # Allow standard numeric types
            # Represent scalar as a tuple expression `(value,)` with key (0,)
            other_expr_str = f"({other!r},)" # Use repr for numbers
            other_keys_tuple = (0,)
        else:
            # Attempt to treat 'other' as a symbolic scalar if it's not a recorder or number
            try:
                # Check if it looks like a symbolic variable or simple expression
                # This is heuristic; ideally, use a dedicated symbolic scalar type
                 if isinstance(other, str) or hasattr(other, '_sympy_'):
                      other_expr_str = f"({other},)" # Assume it can be placed in expr directly
                      other_keys_tuple = (0,)
                 else:
                      return NotImplemented # Cannot handle this type
            except Exception:
                 return NotImplemented

        # Get the operator output (keys_out, func)
        op_output = self._get_operator_output(operator, (other_keys_tuple,))

        # Construct the new expression string
        # Assumes func.__name__ holds the generated function's name
        new_expr = f'{op_output.func.__name__}({self.expr}, {other_expr_str})'

        return TapeRecorder(algebra=self.algebra, expr=new_expr, keys=op_output.keys_out)

    def unary_operator(self, operator: str) -> TapeRecorder:
        """ Apply a unary operator: OP self """
        # Get the operator output (keys_out, func)
        op_output = self._get_operator_output(operator) # No other keys needed

        # Construct the new expression string
        new_expr = f'{op_output.func.__name__}({self.expr})'

        return TapeRecorder(algebra=self.algebra, expr=new_expr, keys=op_output.keys_out)

    # --- Dunder Methods using Operators ---

    # Binary operators
    gp = __mul__ = partialmethod(binary_operator, operator='gp')
    # __rmul__ needs special handling if other is not TapeRecorder
    def __rmul__(self, other: Any) -> TapeRecorder:
        if isinstance(other, (int, float, complex, str)): # Check for types handled by binary_operator
             # Convert scalar 'other' to TapeRecorder and call __mul__
             try:
                  scalar_recorder = TapeRecorder(self.algebra, f"({other!r},)", (0,))
                  return scalar_recorder * self # Calls __mul__
             except TypeError as e:
                  raise TypeError(f"Could not create scalar TapeRecorder for rmul operand: {e}") from e
        return NotImplemented

    sw = __rshift__ = partialmethod(binary_operator, operator='sw')
    # No __rrshift__ needed as order matters

    cp = partialmethod(binary_operator, operator='cp')
    acp = partialmethod(binary_operator, operator='acp')
    ip = __or__ = partialmethod(binary_operator, operator='ip')
    # __ror__ for completeness (IP is often symmetric or defined based on grades)
    def __ror__(self, other: Any) -> TapeRecorder:
        if isinstance(other, (int, float, complex, str)):
             try: scalar_recorder = TapeRecorder(self.algebra, f"({other!r},)", (0,))
             except TypeError as e: raise TypeError(f"Could not create scalar TapeRecorder for ror operand: {e}") from e
             return scalar_recorder | self # Calls __or__
        return NotImplemented

    sp = partialmethod(binary_operator, operator='sp') # Scalar product
    lc = partialmethod(binary_operator, operator='lc') # Left contraction
    rc = partialmethod(binary_operator, operator='rc') # Right contraction

    op = __xor__ = partialmethod(binary_operator, operator='op')
    # __rxor__ for completeness
    def __rxor__(self, other: Any) -> TapeRecorder:
        if isinstance(other, (int, float, complex, str)):
             try: scalar_recorder = TapeRecorder(self.algebra, f"({other!r},)", (0,))
             except TypeError as e: raise TypeError(f"Could not create scalar TapeRecorder for rxor operand: {e}") from e
             return scalar_recorder ^ self # Calls __xor__
        return NotImplemented

    rp = __and__ = partialmethod(binary_operator, operator='rp')
    # __rand__ for completeness
    def __rand__(self, other: Any) -> TapeRecorder:
        if isinstance(other, (int, float, complex, str)):
             try: scalar_recorder = TapeRecorder(self.algebra, f"({other!r},)", (0,))
             except TypeError as e: raise TypeError(f"Could not create scalar TapeRecorder for rand operand: {e}") from e
             return scalar_recorder & self # Calls __and__
        return NotImplemented

    proj = __matmul__ = partialmethod(binary_operator, operator='proj')
    # __rmatmul__ might not make sense depending on proj definition

    add = __add__ = partialmethod(binary_operator, operator='add')
    __radd__ = __add__ # Addition is commutative

    sub = __sub__ = partialmethod(binary_operator, operator='sub')
    def __rsub__(self, other: Any) -> TapeRecorder:
        """ Subtract: other - self """
        if isinstance(other, (int, float, complex, str)):
             try: scalar_recorder = TapeRecorder(self.algebra, f"({other!r},)", (0,))
             except TypeError as e: raise TypeError(f"Could not create scalar TapeRecorder for rsub operand: {e}") from e
             return scalar_recorder - self # Calls __sub__
        return NotImplemented

    div = __truediv__ = partialmethod(binary_operator, operator='div')
    def __rtruediv__(self, other: Any) -> TapeRecorder:
        """ Divide: other / self """
        if isinstance(other, (int, float, complex, str)):
             try: scalar_recorder = TapeRecorder(self.algebra, f"({other!r},)", (0,))
             except TypeError as e: raise TypeError(f"Could not create scalar TapeRecorder for rtruediv operand: {e}") from e
             return scalar_recorder / self # Calls __truediv__
        return NotImplemented


    def __pow__(self, power: int) -> TapeRecorder:
        """ Raise this expression to a non-negative integer power using GP. """
        if not isinstance(power, int):
             raise TypeError(f"Power must be an integer, got {type(power).__name__}")
        if power < 0:
            # Consistent with MultiVector, use inv() for negative powers
            # Requires inv() to be robustly implemented for TapeRecorder
            try:
                 inv_self = self.inv()
                 return inv_self ** (-power)
            except (NotImplementedError, AlgebraError) as e:
                 raise type(e)(f"Cannot raise to negative power {power} as inverse failed: {e}") from e
            except Exception as e: # Catch other potential errors from inv()
                 raise AlgebraError(f"Unexpected error computing inverse for negative power {power}: {e}") from e

        if power == 0:
            # Return a scalar 1 TapeRecorder
            return TapeRecorder(self.algebra, expr='(1,)', keys=(0,))
        if power == 1:
            return self # Return self, no operation needed

        # Exponentiation by squaring using the geometric product (gp)
        res = TapeRecorder(self.algebra, expr='(1,)', keys=(0,)) # Start with identity
        base = self
        p = power
        while p > 0:
             if p % 2 == 1: res = res * base # Use __mul__ (gp)
             base = base * base
             p //= 2
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

    # Methods mirroring MultiVector API
    def dual(self, kind: str = 'auto') -> TapeRecorder:
        """ Compute the dual of this expression. """
        # Mirrors MultiVector.dual logic
        if kind == 'polarity': return self.polarity()
        elif kind == 'hodge': return self.hodge()
        elif kind == 'auto':
            if hasattr(self.algebra, 'r'):
                if self.algebra.r == 0:
                     try: return self.polarity()
                     except (NotImplementedError, AlgebraError): pass # Try hodge if polarity fails/not impl.
                # If r=1 or polarity failed, try hodge
                try: return self.hodge()
                except (NotImplementedError, AlgebraError) as e_hodge:
                     raise AlgebraError(f"Cannot auto-select dual: Hodge failed or not implemented ({e_hodge})") from e_hodge
            raise ValueError("Cannot auto-select dual for this algebra (attribute 'r' missing?).")
        else:
            raise ValueError(f"Unknown dual kind: '{kind}'. Choose 'polarity', 'hodge', or 'auto'.")

    def undual(self, kind: str = 'auto') -> TapeRecorder:
        """ Compute the undual (inverse of dual) of this expression. """
        # Mirrors MultiVector.undual logic
        if kind == 'polarity': return self.unpolarity()
        elif kind == 'hodge': return self.unhodge()
        elif kind == 'auto':
            if hasattr(self.algebra, 'r'):
                if self.algebra.r == 0:
                     try: return self.unpolarity()
                     except (NotImplementedError, AlgebraError): pass
                try: return self.unhodge()
                except (NotImplementedError, AlgebraError) as e_unhodge:
                     raise AlgebraError(f"Cannot auto-select undual: Unhodge failed or not implemented ({e_unhodge})") from e_unhodge
            raise ValueError("Cannot auto-select undual for this algebra (attribute 'r' missing?).")
        else:
            raise ValueError(f"Unknown undual kind: '{kind}'. Choose 'polarity', 'hodge', or 'auto'.")

    def norm(self) -> TapeRecorder:
        """ Compute the norm: sqrt(<self * ~self>_0). """
        # Note: This assumes normsq returns a scalar recorder (key=(0,))
        try:
             normsq_rec = self.normsq() # May raise AlgebraError if x*~x is non-scalar
             # Apply sqrt to the scalar result
             return normsq_rec.sqrt() # sqrt operator should handle scalar recorder
        except (NotImplementedError, AlgebraError) as e:
             raise type(e)(f"Cannot compute norm as normsq/sqrt failed: {e}") from e

    def normalized(self) -> TapeRecorder:
        """ Compute the normalized version: self / self.norm(). """
        try:
             norm_rec = self.norm()
             # Division might fail if norm is symbolically zero or inverse fails
             return self / norm_rec
        except (ZeroDivisionError, NotImplementedError, AlgebraError) as e:
             raise type(e)(f"Cannot normalize TapeRecorder: norm calculation or division failed. {e}") from e

    def compile(self) -> Callable:
        """ Placeholder for compiling the expression string into a function. """
        # (Method unchanged - still NotImplemented)
        raise NotImplementedError("Compilation of TapeRecorder expression not implemented")

    @property
    def e(self) -> TapeRecorder:
        """ Get the scalar component (grade 0) of this expression. """
        # Use grade extraction for consistency
        return self.grade(0)