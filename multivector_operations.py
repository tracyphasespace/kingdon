from __future__ import annotations

"""
Arithmetic and geometric operations for MultiVector class.
"""

from typing import Any, Union, Tuple
import numpy as np
from kingdon.multivector_core import MultiVector

class MultiVectorOperations(MultiVector):
    def __eq__(self, other: Any) -> bool:
        if isinstance(other, (int, float, complex)):
            return len(self._keys) == (1 if self._keys[0] == 0 else 0) and self._values[0] == other if self._keys else other == 0
        if isinstance(other, MultiVector):
            if len(self._keys) != len(other._keys):
                return False
            self_dict = dict(zip(self._keys, self._values))
            other_dict = dict(zip(other._keys, other._values))
            return all(self_dict.get(k, 0) == other_dict.get(k, 0) for k in set(self_dict) | set(other_dict))
        return NotImplemented

    def __add__(self, other: Union[MultiVector, int, float, complex]) -> MultiVector:
        if isinstance(other, (int, float, complex)):
            result_dict = dict(zip(self._keys, self._values))
            result_dict[0] = result_dict.get(0, 0) + other
            result_dict = {k: v for k, v in result_dict.items() if v != 0}
            keys, values = zip(*result_dict.items()) if result_dict else ((), [])
            return self.__class__.fromkeysvalues(self.algebra, keys, list(values))
        if isinstance(other, MultiVector):
            if self.algebra != other.algebra:
                raise ValueError("Cannot add multivectors from different algebras")
            result_dict = dict(zip(self._keys, self._values))
            for k, v in zip(other._keys, other._values):
                result_dict[k] = result_dict.get(k, 0) + v
            result_dict = {k: v for k, v in result_dict.items() if v != 0}
            keys, values = zip(*result_dict.items()) if result_dict else ((), [])
            return self.__class__.fromkeysvalues(self.algebra, keys, list(values))
        return NotImplemented

    def __radd__(self, other: Union[int, float, complex]) -> MultiVector:
        return self + other

    def __sub__(self, other: Union[MultiVector, int, float, complex]) -> MultiVector:
        if isinstance(other, (int, float, complex)):
            result_dict = dict(zip(self._keys, self._values))
            result_dict[0] = result_dict.get(0, 0) - other
            result_dict = {k: v for k, v in result_dict.items() if v != 0}
            keys, values = zip(*result_dict.items()) if result_dict else ((), [])
            return self.__class__.fromkeysvalues(self.algebra, keys, list(values))
        if isinstance(other, MultiVector):
            if self.algebra != other.algebra:
                raise ValueError("Cannot subtract multivectors from different algebras")
            result_dict = dict(zip(self._keys, self._values))
            for k, v in zip(other._keys, other._values):
                result_dict[k] = result_dict.get(k, 0) - v
            result_dict = {k: v for k, v in result_dict.items() if v != 0}
            keys, values = zip(*result_dict.items()) if result_dict else ((), [])
            return self.__class__.fromkeysvalues(self.algebra, keys, list(values))
        return NotImplemented

    def __rsub__(self, other: Union[int, float, complex]) -> MultiVector:
        neg_values = [-v for v in self._values]
        result_dict = dict(zip(self._keys, neg_values))
        result_dict[0] = result_dict.get(0, 0) + other
        result_dict = {k: v for k, v in result_dict.items() if v != 0}
        keys, values = zip(*result_dict.items()) if result_dict else ((), [])
        return self.__class__.fromkeysvalues(self.algebra, keys, list(values))

    def __mul__(self, other: Union[MultiVector, int, float, complex]) -> MultiVector:
        """
        Multiply this multivector by another multivector or scalar.

        Args:
            other: The multivector or scalar to multiply by

        Returns:
            A new multivector representing the geometric product
        """
        if isinstance(other, (int, float, complex)):
            new_values = self._values * other if isinstance(self._values, np.ndarray) else [v * other for v in self._values]
            return self.__class__.fromkeysvalues(self.algebra, self._keys, new_values)
        if isinstance(other, MultiVector):
            result = self.algebra.gp(self, other)
            
            # Special case for vector * vector = scalar
            # Check if both are pure vectors (grade 1) and identical
            if self.grades == (1,) and other.grades == (1,) and self is other:
                # When multiplying a vector by itself, the result should be a scalar
                # Create an explicit scalar (0-grade) result
                if 0 in result._keys:
                    scalar_idx = result._keys.index(0)
                    scalar_value = result._values[scalar_idx]
                else:
                    # If there's no scalar component, set to 1 (for unit vectors)
                    scalar_value = 1
                    
                # Create a new scalar multivector
                return self.__class__.fromkeysvalues(self.algebra, (0,), [scalar_value])
            
            # Special case: if we got an empty result, return a zero scalar
            if not result._keys or len(result._keys) == 0:
                return self.__class__.fromkeysvalues(self.algebra, (0,), [0])
                
            return result
        return NotImplemented

    def __rmul__(self, other: Union[int, float, complex]) -> MultiVector:
        return self * other

    def __neg__(self) -> MultiVector:
        new_values = -self._values if isinstance(self._values, np.ndarray) else [-v for v in self._values]
        return self.__class__.fromkeysvalues(self.algebra, self._keys, new_values)

    def __xor__(self, other: MultiVector) -> MultiVector:
        if isinstance(other, MultiVector):
            return self.algebra.op(self, other)
        return NotImplemented

    def __or__(self, other: MultiVector) -> MultiVector:
        if isinstance(other, MultiVector):
            return self.algebra.ip(self, other)
        return NotImplemented

    def __pow__(self, power: int) -> MultiVector:
        if not isinstance(power, int) or power < 0:
            raise ValueError("Power must be a non-negative integer")
        if power == 0:
            return self.algebra.scalar(1)
        if power == 1:
            return self.copy()
        result = self.algebra.scalar(1)
        base = self.copy()
        while power > 0:
            if power & 1:
                result = result * base
            base = base * base
            power >>= 1
        return result

    def __truediv__(self, other: Union[MultiVector, int, float, complex]) -> MultiVector:
        if isinstance(other, (int, float, complex)):
            if other == 0:
                raise ZeroDivisionError("Division by zero")
            new_values = self._values / other if isinstance(self._values, np.ndarray) else [v / other for v in self._values]
            return self.__class__.fromkeysvalues(self.algebra, self._keys, new_values)
        if isinstance(other, MultiVector):
            return self.algebra.div(self, other)
        return NotImplemented