from __future__ import annotations

"""
Core MultiVector functionality for Geometric Algebra computations.
"""

from typing import Any, Dict, List, Optional, Tuple, Union, TypeVar
import numpy as np
import sympy
from sympy import Expr

try:
    from kingdon.multivector_base import MultiVectorBase
    from kingdon.multivector_indexing import MultiVectorIndexing
    from kingdon.polynomial import RationalPolynomial
except ImportError:
    class MultiVectorBase:
        pass
    class MultiVectorIndexing(MultiVectorBase):
        pass
    class RationalPolynomial:
        pass

T = TypeVar('T', int, float, complex, Expr)

class MultiVector(MultiVectorIndexing):
    """
    Core MultiVector class for Geometric Algebra computations.
    
    Attributes:
        algebra: The associated Geometric Algebra instance
        _values: The coefficients (list or numpy array)
        _keys: The binary keys for basis blades (tuple of integers)
        shape: The shape of the multivector
    """
    
    def __init__(self, algebra: Any, **kwargs: Any) -> None:
        """
        Initialize a MultiVector with named blade coefficients or parameters.
        """
        print(f"MultiVector.__init__: algebra={type(algebra).__name__}, kwargs={kwargs}")
        
        # Initialize core attributes if not already set
        self.algebra = algebra
        
        # Check if _keys and _values are already set (from __new__ or elsewhere)
        if not hasattr(self, '_keys'):
            self._keys = kwargs.get('keys', (0,))
        
        if not hasattr(self, '_values'):
            # Initialize with zeros or provided values
            if 'values' in kwargs:
                values = kwargs['values']
                self._values = values if isinstance(values, (list, np.ndarray)) else [values]
            else:
                self._values = [0] * len(self._keys)
        
        if not hasattr(self, 'shape'):
            # Set shape based on values
            self.shape = self._values.shape if isinstance(self._values, np.ndarray) and self._values.ndim > 1 else ()

    def __getattr__(self, name: str) -> Any:
        if hasattr(self.algebra, 'canon2bin'):
            blade_key = self.algebra.canon2bin.get(name)
            if blade_key is not None:
                return self._values[self._keys.index(blade_key)] if blade_key in self._keys else 0
        raise AttributeError(f"'{type(self).__name__}' object has no attribute '{name}'")

    def __len__(self) -> int:
        return len(self._values)

    def __str__(self) -> str:
        if not self._values or len(self._values) == 0:
            return "0"
        def format_val(val: Any) -> str:
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
                components.append(f"{format_val(val)}{blade}" if val not in (1, -1) else f"-{blade}" if val == -1 else blade)
        return " + ".join(components).replace("+ -", "- ") if components else "0"

    def __repr__(self) -> str:
        return f"MultiVector({self.algebra}, {dict(self.items())})"
        
    def keys(self) -> Tuple[int, ...]:
        """
        Return the binary keys of the multivector.
        
        Each key is an integer representing a basis blade in the algebra.
        
        Returns:
            Tuple[int, ...]: Binary keys corresponding to basis blades
        """
        return self._keys

    def values(self) -> Union[List[Any], np.ndarray]:
        """
        Return the list or array of coefficients.
        
        Returns:
            Union[List[Any], np.ndarray]: Coefficients for each basis blade
        """
        return self._values

    def items(self) -> List[Tuple[int, Any]]:
        """
        Return a list of (key, value) pairs.
        
        Returns:
            List[Tuple[int, Any]]: Pairs of (key, coefficient)
        """
        return list(zip(self._keys, self._values))