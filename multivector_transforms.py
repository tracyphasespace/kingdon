from __future__ import annotations

"""
Transformation methods for MultiVector class.
"""

import math
from typing import Union
import numpy as np
from kingdon.multivector_operations import MultiVectorOperations

class MultiVectorTransforms(MultiVectorOperations):
    def copy(self) -> MultiVectorTransforms:
        values_copy = self._values.copy() if isinstance(self._values, np.ndarray) else list(self._values)
        return self.__class__.fromkeysvalues(self.algebra, self._keys, values_copy)

    def project_grade(self, grade: int) -> MultiVectorTransforms:
        matching_keys = [k for k in self._keys if bin(k).count('1') == grade]
        if not matching_keys:
            return self.__class__.fromkeysvalues(self.algebra, (), [])
        matching_indices = [self._keys.index(k) for k in matching_keys]
        matching_values = self._values[matching_indices] if isinstance(self._values, np.ndarray) else [self._values[i] for i in matching_indices]
        return self.__class__.fromkeysvalues(self.algebra, tuple(matching_keys), matching_values)

    def grade(self, *grades: int) -> MultiVectorTransforms:
        if not grades:
            return self.__class__.fromkeysvalues(self.algebra, (), [])
        matching_keys = [k for k in self._keys if bin(k).count('1') in grades]
        if not matching_keys:
            return self.__class__.fromkeysvalues(self.algebra, (), [])
        matching_indices = [self._keys.index(k) for k in matching_keys]
        matching_values = self._values[matching_indices] if isinstance(self._values, np.ndarray) else [self._values[i] for i in matching_indices]
        return self.__class__.fromkeysvalues(self.algebra, tuple(matching_keys), matching_values)

    def scalar(self) -> MultiVectorTransforms:
        return self.project_grade(0)

    def vector(self) -> MultiVectorTransforms:
        return self.project_grade(1)

    def bivector(self) -> MultiVectorTransforms:
        return self.project_grade(2)

    def reverse(self) -> MultiVectorTransforms:
        # Implement reverse directly to avoid relying on algebra.reverse
        result = {}
        for k, v in self.items():
            grade = bin(k).count('1')  # Calculate the grade
            sign = -1 if grade % 4 in (2, 3) else 1  # Grades 2 and 3 (mod 4) change sign
            result[k] = -v if sign < 0 else v  # Apply sign change
            
        # Convert result to a proper multivector
        keys, values = zip(*result.items()) if result else ((), [])
        return self.__class__.fromkeysvalues(self.algebra, keys, list(values))

    def involute(self) -> MultiVectorTransforms:
        return self.algebra.involute(self)

    def conjugate(self) -> MultiVectorTransforms:
        return self.algebra.conjugate(self)

    def __invert__(self) -> MultiVectorTransforms:
        return self.reverse()

    def dual(self) -> MultiVectorTransforms:
        return self.algebra.dual(self)

    def inverse(self) -> MultiVectorTransforms:
        return self.algebra.inv(self)

    def norm(self) -> float:
        """
        Compute the norm (magnitude) of this multivector.
        
        The norm is computed as sqrt(|a * ~a|) where ~a is the reverse.
        For vectors this is the usual Euclidean norm.
        
        Returns:
            The norm as a float value
            
        Raises:
            ZeroDivisionError: If the multivector has zero norm
        """
        # Calculate |A * ~A| = |A|²
        norm_sq = (self * self.reverse()).scalar()
        if not norm_sq:
            return 0.0
            
        # Extract the scalar value
        if hasattr(norm_sq, 'e'):
            norm_sq_val = norm_sq.e
        else:
            # If e attribute is not available, try to access the first value
            norm_sq_val = norm_sq._values[0] if len(norm_sq._values) > 0 else 0
            
        # Take the square root
        if norm_sq_val < 0:
            # Handle negative squared norm (can happen with certain signatures)
            return math.sqrt(abs(norm_sq_val))
        else:
            return math.sqrt(norm_sq_val)

    def normalized(self) -> MultiVectorTransforms:
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

    def __array__(self) -> np.ndarray:
        return self._values if isinstance(self._values, np.ndarray) else np.array(self._values)

__all__ = ['MultiVectorTransforms']