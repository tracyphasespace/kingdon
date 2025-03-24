from __future__ import annotations

import numpy as np
import sympy
from dataclasses import dataclass, field
from functools import cached_property
from typing import (
    Any, 
    Callable, 
    Dict, 
    Generator, 
    List, 
    Optional, 
    Tuple, 
    Union, 
    Set, 
    TypeVar
)

# Try to import RationalPolynomial without causing circular imports
try:
    from kingdon.polynomial import RationalPolynomial
except ImportError:
    # Create a placeholder class if it can't be imported
    class RationalPolynomial:
        pass

# Type variable for scalar values
T = TypeVar('T', int, float, complex, sympy.Expr, RationalPolynomial)

@dataclass(init=False)
class MultiVectorBase:
    """
    Base class for MultiVector representation in Geometric Algebra.
    
    Provides core functionality for storing and manipulating multivector components.
    
    Attributes:
        algebra: The associated Geometric Algebra instance
        _values: Coefficients for the multivector components (list or numpy array)
        _keys: Binary keys corresponding to basis blades (tuple of integers)
        shape: Shape of the multivector (empty tuple for scalar, specific shape for array-valued)
    """
    algebra: Any
    _values: Union[List[Any], np.ndarray] = field(default_factory=list)
    _keys: Tuple[int, ...] = field(default_factory=tuple)
    shape: Tuple[int, ...] = field(default_factory=tuple)

    @property
    def d(self) -> int:
        """
        Get the dimension of the algebra.
        
        Returns:
            int: Total number of dimensions in the algebra
        """
        return self.algebra.d if hasattr(self.algebra, 'd') else 0

    @property
    def coeffs(self) -> Dict[int, Any]:
        """
        Return a dictionary mapping basis blade keys to their coefficients.
        
        For array-valued multivectors, this returns the coefficients for the first element.
        
        Returns:
            Dict[int, Any]: Dictionary mapping binary keys to their coefficients
        """
        if isinstance(self._values, np.ndarray) and self._values.ndim > 1:
            # For array-valued multivectors, return the first element's coefficients
            return dict(zip(self._keys, self._values.flat[:len(self._keys)]))
        return dict(zip(self._keys, self._values))

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

    def items(self) -> Generator[Tuple[int, Any], None, None]:
        """
        Return an iterator of (key, value) pairs.
        
        Yields:
            Tuple[int, Any]: Pairs of (key, coefficient)
        """
        yield from zip(self._keys, self._values)

    def __len__(self) -> int:
        """
        Return the number of components in the multivector.
        
        Returns:
            int: Number of (key, coefficient) pairs
        """
        return len(self._values)

    @cached_property
    def type_number(self) -> int:
        """
        Compute a unique type number based on present basis blades.
        
        Returns:
            int: The type number as an integer.
        """
        bitstr = "".join(
            "1" if i in self.keys() else "0" 
            for i in reversed(self.algebra.canon2bin.values())
        )
        return int(bitstr, 2)

    @cached_property
    def grades(self) -> Tuple[int, ...]:
        """
        Return the sorted tuple of grades present in the multivector.
        
        Returns:
            Tuple[int, ...]: Grades present in the multivector.
        """
        if not self.keys():
            return tuple()
            
        # Filter out non-zero coefficients
        grades_dict = {}
        for k, v in zip(self._keys, self._values):
            grade = bin(k).count("1")
            # Only include non-zero components
            is_zero = (
                (isinstance(v, (int, float, complex)) and v == 0) or
                (isinstance(v, np.ndarray) and np.all(v == 0))
            )
            if not is_zero:
                grades_dict[grade] = True
                
        # If no non-zero coefficients found, ensure we still return grade 0 for zero scalar
        if not grades_dict and 0 in [bin(k).count("1") for k in self.keys()]:
            grades_dict[0] = True
            
        return tuple(sorted(grades_dict.keys()))

    @cached_property
    def issymbolic(self) -> bool:
        """
        Determine if the multivector contains symbolic coefficients.
        
        Returns:
            bool: True if any coefficient is symbolic, False otherwise.
        """
        symbol_classes = (sympy.Expr, RationalPolynomial)
        if hasattr(self.algebra, 'codegen_symbolcls') and self.algebra.codegen_symbolcls:
            symbolcls = self.algebra.codegen_symbolcls
            symbol_classes = (*symbol_classes, 
                symbolcls.__self__ if hasattr(symbolcls, "__self__") else symbolcls)
            
        return any(isinstance(v, symbol_classes) for v in self.values())

    @cached_property
    def free_symbols(self) -> Set:
        """
        Return the set of free symbols present in the multivector.
        
        Returns:
            Set: A set of sympy Symbols.
        """
        symbols = set()
        for v in self.values():
            if hasattr(v, "free_symbols"):
                symbols |= v.free_symbols
        return symbols