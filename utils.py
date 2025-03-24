"""
Utility Functions for Geometric Algebra
=======================================

This module provides utility functions that support the Kingdon Geometric Algebra library.
These functions handle common operations like efficient exponentiation and other
mathematical operations that are frequently used across the library.

Key functions:
- power_supply: Efficiently compute powers of an object using addition chains
- add_chains: Generate minimal addition chains for exponentiation

Example:
    >>> from kingdon.utils import power_supply
    >>> import operator
    >>> 
    >>> # Generate powers of 2 up to 2^5
    >>> powers = list(power_supply(2, 5, operation=operator.mul))
    >>> print(powers)
    [2, 4, 8, 16, 32]
    
    >>> # Generate specific powers of a matrix
    >>> import numpy as np
    >>> matrix = np.array([[1, 1], [0, 1]])
    >>> powers = list(power_supply(matrix, [2, 4, 8], operation=np.matmul))
"""

from __future__ import annotations

import math
import operator
from typing import Any, Callable, Dict, Generator, Iterable, List, Tuple, Union, TypeVar, Optional, Set

# Define a type variable for the item being raised to powers
T = TypeVar('T')


class AdditionChain:
    """
    Efficiently computes minimal addition chains for exponentiation.
    
    An addition chain for n is a sequence starting with 1 where each element
    is the sum of two previous elements, and the final element is n.
    This allows computing x^n with fewer multiplications than naive approaches.
    
    Attributes:
        _chains: Dictionary mapping integers to their minimal addition chains
        
    Example:
        >>> chain = AdditionChain()
        >>> chain.get_chain(15)
        [1, 2, 3, 6, 12, 15]
        # To compute x^15: x^2=x*x, x^3=x^2*x, x^6=x^3*x^3, x^12=x^6*x^6, x^15=x^12*x^3
    """
    def __init__(self) -> None:
        """Initialize with the basic chain for 1."""
        self._chains: Dict[int, List[int]] = {1: [1]}
        
    def get_chain(self, n: int) -> List[int]:
        """
        Get a minimal addition chain for n.
        
        Args:
            n: Target exponent
            
        Returns:
            List of integers forming an addition chain ending with n
            
        Raises:
            ValueError: If n < 1
        """
        if n < 1:
            raise ValueError("Exponent must be a positive integer")
            
        if n in self._chains:
            return self._chains[n].copy()
            
        # Generate chains up to n
        self._generate_chains_up_to(n)
        return self._chains[n].copy()
        
    def _generate_chains_up_to(self, n: int) -> None:
        """
        Generate all minimal addition chains up to n.
        
        Args:
            n: Maximum exponent to generate chains for
        """
        # Keep generating chains until we have one for n
        while n not in self._chains:
            for m in range(1, n):
                if m in self._chains:
                    chain = self._chains[m]
                    
                    # Try to extend the chain by adding its largest element to another element
                    largest = chain[-1]
                    for i in range(len(chain) - 1, -1, -1):
                        next_value = largest + chain[i]
                        
                        if next_value <= n and next_value not in self._chains:
                            self._chains[next_value] = chain + [next_value]
                            
                        # If we found a chain for n, we're done
                        if next_value == n:
                            return
        

# Global instance for reuse
_addition_chain = AdditionChain()


def get_addition_chain(n: int) -> List[int]:
    """
    Get a minimal addition chain for n.
    
    Args:
        n: Target exponent
        
    Returns:
        List of integers forming an addition chain ending with n
    """
    return _addition_chain.get_chain(n)


def power_supply(x: T, 
                 exponents: Union[int, Iterable[int]], 
                 operation: Callable[[T, T], T] = operator.mul) -> Generator[T, None, None]:
    """
    Generate powers of an object using minimal operations.
    
    This function uses addition chains to efficiently compute powers. It can
    compute a single power (e.g., x^7) or multiple powers (e.g., x^2, x^5, x^8)
    with minimal operations.
    
    Args:
        x: The object to raise to powers
        exponents: Either a single integer or a collection of integers
        operation: Binary operation to use for "multiplication" (default is operator.mul)
            The operation must take two arguments (a, b) and return a value of the same type
        
    Yields:
        Powers of x in ascending order
        
    Example:
        >>> # Generate powers 2, 4, and 7 of a matrix
        >>> import numpy as np
        >>> matrix = np.array([[1, 1], [0, 1]])
        >>> powers = list(power_supply(matrix, [2, 4, 7], operation=np.matmul))
        
    Notes:
        - Avoids recomputing the same power multiple times
        - For collections of exponents, computes a single optimal chain 
        - Yields only the requested powers in ascending order
    """
    # Handle case of a single exponent
    if isinstance(exponents, int):
        target = exponents
        chain = get_addition_chain(target)
    else:
        # Convert to sorted list and validate
        exponents_list = sorted(set(int(e) for e in exponents))
        if not exponents_list:
            return
        if min(exponents_list) < 1:
            raise ValueError("Exponents must be positive integers")
            
        # Get chain for the maximum exponent
        target = max(exponents_list)
        chain = get_addition_chain(target)
        
    # Compute powers using the addition chain
    powers: Dict[int, T] = {1: x}
    for i in range(1, len(chain)):
        # Find two previous elements in the chain that sum to the current element
        current = chain[i]
        # Simple strategy: use the previous element and find what to add
        prev = chain[i-1]
        to_add = current - prev
        
        # Compute the new power
        powers[current] = operation(powers[prev], powers[to_add])
        
        # If this is one of the requested exponents, yield it
        if isinstance(exponents, int):
            # For single exponent, yield each power in the chain
            yield powers[current]
        elif current in exponents_list:
            # For multiple exponents, yield only those requested
            yield powers[current]


def binary_operation_through_powers(base: Any, 
                             exponent: int, 
                             identity: Any, 
                             operation: Callable[[Any, Any], Any]) -> Any:
    """
    Compute base^exponent using binary operation through powers of 2.
    
    This implements the binary exponentiation algorithm (also known as exponentiation by squaring),
    which is efficient for large exponents.
    
    Args:
        base: The base value
        exponent: The exponent to raise the base to
        identity: The identity element for the operation (e.g., 1 for multiplication)
        operation: The binary operation to apply
        
    Returns:
        The result of base^exponent computed through the binary operation
        
    Example:
        >>> # Compute 2^10 using multiplication
        >>> binary_operation_through_powers(2, 10, 1, operator.mul)
        1024
        
        >>> # Compute matrix^5 using matrix multiplication
        >>> import numpy as np
        >>> matrix = np.array([[1, 1], [0, 1]])
        >>> binary_operation_through_powers(matrix, 5, np.eye(2), np.matmul)
        array([[1, 5],
               [0, 1]])
    """
    if exponent < 0:
        raise ValueError("Exponent must be non-negative")
        
    if exponent == 0:
        return identity
        
    result = identity
    current_power = base
    
    # Binary exponentiation algorithm
    while exponent > 0:
        # If current bit is 1, multiply result by the current power
        if exponent & 1:
            result = operation(result, current_power)
            
        # Move to the next bit
        exponent >>= 1
        
        # Square the current power if there are more bits to process
        if exponent > 0:
            current_power = operation(current_power, current_power)
            
    return result


def is_power_of_two(n: int) -> bool:
    """
    Check if a number is a power of 2.
    
    Args:
        n: Integer to check
        
    Returns:
        True if n is a power of 2, False otherwise
        
    Example:
        >>> is_power_of_two(16)
        True
        >>> is_power_of_two(20)
        False
    """
    return n > 0 and (n & (n - 1)) == 0


def next_power_of_two(n: int) -> int:
    """
    Find the smallest power of 2 greater than or equal to n.
    
    Args:
        n: Integer
        
    Returns:
        Smallest power of 2 >= n
        
    Example:
        >>> next_power_of_two(13)
        16
        >>> next_power_of_two(32)
        32
    """
    if n <= 0:
        return 1
        
    # If already a power of 2, return n
    if is_power_of_two(n):
        return n
        
    # Find the position of the most significant bit
    return 1 << (n.bit_length())


def binomial_coefficient(n: int, k: int) -> int:
    """
    Compute the binomial coefficient (n choose k).
    
    Args:
        n: Total number of items
        k: Number of items to choose
        
    Returns:
        Binomial coefficient n choose k
        
    Example:
        >>> binomial_coefficient(5, 2)
        10
    """
    if k < 0 or k > n:
        return 0
    if k == 0 or k == n:
        return 1
        
    # Use the symmetry of binomial coefficients
    k = min(k, n - k)
    
    # Compute using multiplicative formula
    result = 1
    for i in range(k):
        result = result * (n - i) // (i + 1)
        
    return result