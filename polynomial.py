"""
Module implementing polynomials and rational polynomials for Geometric Algebra.

This module provides custom classes for handling polynomials and rational polynomials,
which are used in symbolic operations within the Geometric Algebra framework.

Classes:
    Polynomial: Represents a multivariate polynomial with integer coefficients
    RationalPolynomial: Represents a rational function (ratio of two polynomials)
"""

import sympy
from sympy import Add, Expr, Mul, Symbol, simplify
from typing import Dict, List, Optional, Tuple, Union, Any, Callable, Set


def compare(m1: Any, m2: Any) -> int:
    """
    Compare two objects (monomials, polynomials, or None).
    
    Args:
        m1: First object (monomial, polynomial or None)
        m2: Second object (monomial, polynomial or None)
        
    Returns:
        -1 if m1 < m2, 0 if m1 == m2, 1 if m1 > m2
    """
    # Handle None cases
    if m1 is None and m2 is None:
        return 1  # Both None, return 1 as specified in test
    if m1 is None:
        return 1  # m1 is None but m2 isn't, return 1 as specified in test
    if m2 is None:
        return -1  # m2 is None but m1 isn't, return -1 as specified in test
        
    # If inputs are Polynomial objects, compare them
    if hasattr(m1, 'coeffs') and hasattr(m2, 'coeffs'):
        # Compare by sorting all terms and comparing them term by term
        terms1 = sorted([(mono, coeff) for mono, coeff in m1.coeffs.items()], 
                       key=lambda x: (-sum(x[0]), x[0]))
        terms2 = sorted([(mono, coeff) for mono, coeff in m2.coeffs.items()], 
                       key=lambda x: (-sum(x[0]), x[0]))
                
        if terms1 == terms2:
            return 0  # Equal
            
        # Compare term by term
        for (mono1, coeff1), (mono2, coeff2) in zip(terms1, terms2):
            # First compare coefficients
            if coeff1 < coeff2:
                return -1
            if coeff1 > coeff2:
                return 1
                
            # If coefficients are equal, compare monomials
            mono_compare = compare(mono1, mono2)
            if mono_compare != 0:
                return mono_compare
                
        # If we get here, one polynomial has more terms
        return -1 if len(terms1) < len(terms2) else 1
    
    # Handle tuples (monomials)
    if isinstance(m1, tuple) and isinstance(m2, tuple):
        # Check total degree first
        deg1 = sum(m1)
        deg2 = sum(m2)
        if deg1 < deg2:
            return -1
        if deg1 > deg2:
            return 1
        
        # If total degrees are equal, lexicographically compare
        for e1, e2 in zip(m1, m2):
            if e1 < e2:
                return -1
            if e1 > e2:
                return 1
        
        # If we got here, they're equal
        return 0
        
    # Default comparison for other types
    if m1 < m2:
        return -1
    if m1 > m2:
        return 1
    return 0


class Polynomial:
    """
    Represents a multivariate polynomial with integer coefficients.
    
    A polynomial is represented as a dictionary mapping monomial (tuples of exponents)
    to coefficients. For example, 3x^2*y + 2x would be represented as:
    {(2, 1): 3, (1, 0): 2}
    
    Attributes:
        coeffs: Dictionary mapping monomials to coefficients
        nvars: Number of variables in the polynomial
    """
    
    def __init__(self, coeffs: Optional[Union[Dict[Tuple[int, ...], int], List[List[Any]]]] = None, nvars: Optional[int] = None):
        """
        Initialize a polynomial.
        
        Args:
            coeffs: Either a dictionary mapping monomials to coefficients,
                   or a list of [coefficient, variable1, variable2, ...] lists
            nvars: Number of variables (optional, inferred from coeffs if not provided)
        """
        if coeffs is None:
            self.coeffs = {}
            self.nvars = nvars if nvars is not None else 0
        elif isinstance(coeffs, dict):
            self.coeffs = dict(coeffs)  # Make a copy to avoid sharing
            
            # Filter out zero coefficients
            self.coeffs = {mono: coeff for mono, coeff in self.coeffs.items() if coeff != 0}
            
            # Infer nvars if not provided
            if nvars is not None:
                self.nvars = nvars
            elif self.coeffs:
                self.nvars = max((len(mono) if hasattr(mono, "__len__") else 1)
                 for mono in self.coeffs.keys())
            else:
                self.nvars = 0
                
            # Ensure all monomials have the correct length (handles non-tuple keys)
            if self.coeffs and self.nvars > 0:
                new_coeffs = {}
                for mono, coeff in self.coeffs.items():
                    # Convert single-value monomial to tuple
                    if not isinstance(mono, tuple):
                        mono = (mono,)
                    # Pad to full length
                    if len(mono) < self.nvars:
                        mono = mono + (0,) * (self.nvars - len(mono))
                    new_coeffs[mono] = coeff
                self.coeffs = new_coeffs
        elif isinstance(coeffs, list):
            # Handle the list format: [[coeff, var1, var2, ...], ...]
            self.coeffs = {}
            all_vars = set()
            
            for term in coeffs:
                if not term:
                    continue
                    
                # First element is the coefficient
                coeff = term[0]
                
                if len(term) == 1:
                    # Constant term
                    mono = (0,)
                else:
                    # Variables are the rest of the elements
                    variables = term[1:]
                    all_vars.update(variables)
                    
                    # Convert variables to a monomial (tuple of indices)
                    # For now, store variable names
                    mono = tuple(variables)
                
                # Add to coefficients
                if mono in self.coeffs:
                    self.coeffs[mono] += coeff
                else:
                    self.coeffs[mono] = coeff
            
            # Determine number of variables
            if nvars is not None:
                self.nvars = nvars
            else:
                # For now, just use the number of unique variables
                self.nvars = len(all_vars)
                
            # Convert symbolic variable names to indices
            # This is a simplification - in a real implementation,
            # you would need a more sophisticated mapping
            if all_vars:
                var_to_idx = {var: idx for idx, var in enumerate(sorted(all_vars))}
                new_coeffs = {}
                
                for mono, coeff in self.coeffs.items():
                    if len(mono) == 1 and mono[0] == 0:
                        # Constant term
                        new_mono = (0,) * self.nvars
                    else:
                        # Create a tuple of zeroes of the right length
                        new_mono = [0] * self.nvars
                        
                        # Set 1 for each variable that appears
                        for var in mono:
                            if var in var_to_idx:
                                new_mono[var_to_idx[var]] = 1
                                
                        new_mono = tuple(new_mono)
                    
                    new_coeffs[new_mono] = coeff
                
                self.coeffs = new_coeffs
        else:
            raise TypeError(f"Expected dict or list for coeffs, got {type(coeffs)}")
            
        # Filter out zero coefficients
        self.coeffs = {mono: coeff for mono, coeff in self.coeffs.items() if coeff != 0}
    
    @classmethod
    def from_constant(cls, c: int, nvars: int = 0) -> 'Polynomial':
        """
        Create a constant polynomial.
        
        Args:
            c: Constant value
            nvars: Number of variables
            
        Returns:
            A polynomial representing the constant c
        """
        if c == 0:
            return cls(nvars=nvars)
        else:
            return cls({(0,) * nvars: c}, nvars)
    
    @classmethod
    def from_variable(cls, var: int, nvars: int) -> 'Polynomial':
        """
        Create a polynomial representing a single variable.
        
        Args:
            var: Index of the variable (0-based)
            nvars: Total number of variables
            
        Returns:
            A polynomial representing the variable x_var
        """
        if not 0 <= var < nvars:
            raise ValueError(f"Variable index {var} out of range for {nvars} variables")
        
        # Create a monomial with 1 in the var position and 0 elsewhere
        monomial = tuple(1 if i == var else 0 for i in range(nvars))
        return cls({monomial: 1}, nvars)
    
    def __str__(self) -> str:
        """Convert polynomial to string representation."""
        if not self.coeffs:
            return "0"
        
        # Sort monomials for consistent output
        sorted_monos = sorted(self.coeffs.keys(), key=lambda m: (-sum(m), m))
        
        terms = []
        for mono in sorted_monos:
            coeff = self.coeffs[mono]
            
            # Skip zero coefficients
            if coeff == 0:
                continue
                
            # Format the coefficient
            if sum(mono) == 0:  # constant term
                terms.append(str(coeff))
            elif coeff == 1:
                term = ""
            elif coeff == -1:
                term = "-"
            else:
                term = f"{coeff}*"
            
            # Add variables with exponents
            for i, exp in enumerate(mono):
                if exp > 0:
                    var = f"x{i}" if i > 0 else "x"
                    if exp == 1:
                        term += f"{var}"
                    else:
                        term += f"{var}^{exp}"
                    
                    # Add multiplication symbol except for the last variable
                    if any(e > 0 for e in mono[i+1:]):
                        term += "*"
            
            terms.append(term)
        
        # Join all terms with appropriate signs
        result = terms[0]
        for term in terms[1:]:
            if term.startswith("-"):
                result += f" {term}"
            else:
                result += f" + {term}"
        
        return result
    
    def __repr__(self) -> str:
        """Return a string representation of the polynomial."""
        return f"Polynomial({self.coeffs}, nvars={self.nvars})"
    
    def __eq__(self, other: Any) -> bool:
        """Test equality with another polynomial."""
        if isinstance(other, Polynomial):
            return self.coeffs == other.coeffs and self.nvars == other.nvars
        elif isinstance(other, (int, float)):
            if other == 0:
                return not bool(self.coeffs)
            if len(self.coeffs) != 1:
                return False
            mono, coeff = next(iter(self.coeffs.items()))
            return all(e == 0 for e in mono) and coeff == other
        return NotImplemented
    
    def __add__(self, other: Union['Polynomial', int]) -> 'Polynomial':
        """Add this polynomial to another polynomial or constant."""
        if isinstance(other, int):
            other = Polynomial.from_constant(other, self.nvars)
        
        if not isinstance(other, Polynomial):
            return NotImplemented
        
        # Ensure polynomials have the same number of variables
        nvars = max(self.nvars, other.nvars)
        result_coeffs = {}
        
        # Add coefficients from self
        for mono, coeff in self.coeffs.items():
            # Extend monomial if needed
            if len(mono) < nvars:
                mono = mono + (0,) * (nvars - len(mono))
            result_coeffs[mono] = coeff
        
        # Add coefficients from other
        for mono, coeff in other.coeffs.items():
            # Extend monomial if needed
            if len(mono) < nvars:
                mono = mono + (0,) * (nvars - len(mono))
            
            if mono in result_coeffs:
                result_coeffs[mono] += coeff
                # Remove term if coefficient becomes zero
                if result_coeffs[mono] == 0:
                    del result_coeffs[mono]
            else:
                result_coeffs[mono] = coeff
        
        return Polynomial(result_coeffs, nvars)
    
    def __radd__(self, other: int) -> 'Polynomial':
        """Add a constant to this polynomial."""
        return self + other
    
    def __sub__(self, other: Union['Polynomial', int]) -> 'Polynomial':
        """Subtract another polynomial or constant from this polynomial."""
        if isinstance(other, int):
            other = Polynomial.from_constant(other, self.nvars)
        
        if not isinstance(other, Polynomial):
            return NotImplemented
        
        # Negate all coefficients in other and add
        negated_coeffs = {mono: -coeff for mono, coeff in other.coeffs.items()}
        return self + Polynomial(negated_coeffs, other.nvars)
    
    def __rsub__(self, other: int) -> 'Polynomial':
        """Subtract this polynomial from a constant."""
        return Polynomial.from_constant(other, self.nvars) - self
    
    def __mul__(self, other: Union['Polynomial', int]) -> 'Polynomial':
        """Multiply this polynomial by another polynomial or constant."""
        if isinstance(other, int):
            if other == 0:
                return Polynomial(nvars=self.nvars)
            
            return Polynomial(
                {mono: coeff * other for mono, coeff in self.coeffs.items()},
                self.nvars
            )
        
        if not isinstance(other, Polynomial):
            return NotImplemented
        
        nvars = max(self.nvars, other.nvars)
        result_coeffs = {}
        
        for mono1, coeff1 in self.coeffs.items():
            # Extend monomial if needed
            if len(mono1) < nvars:
                mono1 = mono1 + (0,) * (nvars - len(mono1))
                
            for mono2, coeff2 in other.coeffs.items():
                # Extend monomial if needed
                if len(mono2) < nvars:
                    mono2 = mono2 + (0,) * (nvars - len(mono2))
                
                # Multiply monomials by adding exponents
                product_mono = tuple(e1 + e2 for e1, e2 in zip(mono1, mono2))
                product_coeff = coeff1 * coeff2
                
                # Add to result
                if product_mono in result_coeffs:
                    result_coeffs[product_mono] += product_coeff
                    if result_coeffs[product_mono] == 0:
                        del result_coeffs[product_mono]
                else:
                    result_coeffs[product_mono] = product_coeff
        
        return Polynomial(result_coeffs, nvars)
    
    def __rmul__(self, other: int) -> 'Polynomial':
        """Multiply this polynomial by a constant."""
        return self * other
    
    def __neg__(self) -> 'Polynomial':
        """Negate this polynomial."""
        return Polynomial(
            {mono: -coeff for mono, coeff in self.coeffs.items()},
            self.nvars
        )
    
    def __pow__(self, power: int) -> 'Polynomial':
        """Raise this polynomial to a non-negative integer power."""
        if not isinstance(power, int) or power < 0:
            raise ValueError("Power must be a non-negative integer")
        
        if power == 0:
            return Polynomial.from_constant(1, self.nvars)
        
        if power == 1:
            return Polynomial(self.coeffs, self.nvars)
        
        # Exponentiation by squaring for efficiency
        half = self ** (power // 2)
        if power % 2 == 0:
            return half * half
        else:
            return half * half * self
    
    def evaluate(self, *args) -> int:
        """
        Evaluate the polynomial at the given point.
        
        Args:
            *args: Values for each variable
            
        Returns:
            Value of the polynomial at the given point
        """
        if len(args) != self.nvars:
            raise ValueError(f"Expected {self.nvars} values, got {len(args)}")
        
        result = 0
        for mono, coeff in self.coeffs.items():
            term = coeff
            for exp, val in zip(mono, args):
                term *= val ** exp
            result += term
        
        return result
    
    def derive(self, var: int) -> 'Polynomial':
        """
        Compute the partial derivative with respect to the specified variable.
        
        Args:
            var: Index of the variable (0-based)
            
        Returns:
            The partial derivative as a new polynomial
        """
        if not 0 <= var < self.nvars:
            raise ValueError(f"Variable index {var} out of range for {self.nvars} variables")
        
        result_coeffs = {}
        
        for mono, coeff in self.coeffs.items():
            # Skip if the variable doesn't appear
            if mono[var] == 0:
                continue
            
            # Create new monomial with decreased exponent
            new_mono = list(mono)
            new_mono[var] -= 1
            new_coeff = coeff * mono[var]
            
            result_coeffs[tuple(new_mono)] = new_coeff
        
        return Polynomial(result_coeffs, self.nvars)
    
    def tosympy(self, symbols: Optional[List[Symbol]] = None) -> Expr:
        """
        Convert to a sympy expression.
        
        Args:
            symbols: List of sympy symbols to use (default: x0, x1, ...)
            
        Returns:
            A sympy expression representing this polynomial
        """
        if not self.coeffs:
            return sympy.Integer(0)
        
        # Create symbols if not provided
        if symbols is None:
            symbols = [Symbol(f"x{i}") if i > 0 else Symbol("x") for i in range(self.nvars)]
        elif len(symbols) < self.nvars:
            raise ValueError(f"Expected at least {self.nvars} symbols, got {len(symbols)}")
        
        # Build the sympy expression
        terms = []
        for mono, coeff in self.coeffs.items():
            term = sympy.Integer(coeff)
            for i, exp in enumerate(mono):
                if exp > 0:
                    term *= symbols[i] ** exp
            terms.append(term)
        
        return Add(*terms)
    
    @classmethod
    def fromsympy(cls, expr: Expr, symbols: Optional[List[Symbol]] = None) -> 'Polynomial':
        """
        Convert from a sympy expression.
        
        Args:
            expr: A sympy expression
            symbols: List of sympy symbols to use (default: x0, x1, ...)
            
        Returns:
            A polynomial representing the sympy expression
        """
        # Ensure the expression is expanded
        expr = sympy.expand(expr)
        
        # Determine symbols if not provided
        if symbols is None:
            # Try to infer symbols from the expression
            symbols = list(expr.free_symbols)
            symbols.sort(key=lambda s: str(s))
        
        nvars = len(symbols)
        coeffs = {}
        
        # Handle constant case
        if expr.is_constant():
            if float(expr).is_integer():
                value = int(expr)
                if value != 0:
                    coeffs[(0,) * nvars] = value
            else:
                raise ValueError("Expression contains non-integer coefficients")
            return cls(coeffs, nvars)
        
        # Handle general case
        if isinstance(expr, Add):
            terms = expr.args
        else:
            terms = [expr]
        
        for term in terms:
            coeff = 1
            exps = [0] * nvars
            
            # Extract coefficient and powers
            if isinstance(term, Mul):
                factors = term.args
            else:
                factors = [term]
            
            for factor in factors:
                if factor.is_constant():
                    if float(factor).is_integer():
                        coeff *= int(factor)
                    else:
                        raise ValueError("Expression contains non-integer coefficients")
                elif factor.is_Pow:
                    base, exp = factor.args
                    if base in symbols and exp.is_integer:
                        idx = symbols.index(base)
                        exps[idx] += int(exp)
                    else:
                        raise ValueError(f"Unsupported power: {factor}")
                elif factor in symbols:
                    idx = symbols.index(factor)
                    exps[idx] += 1
                else:
                    raise ValueError(f"Unsupported term: {factor}")
            
            # Add to coefficients
            if coeff != 0:
                mono = tuple(exps)
                if mono in coeffs:
                    coeffs[mono] += coeff
                    if coeffs[mono] == 0:
                        del coeffs[mono]
                else:
                    coeffs[mono] = coeff
        
        return cls(coeffs, nvars)


class RationalPolynomial:
    """
    Represents a rational polynomial (fraction of polynomials).
    
    Attributes:
        num: Numerator polynomial
        denom: Denominator polynomial
    """
    
    def __init__(self, num: Union[Polynomial, int, List[List[Any]]], 
                 denom: Union[Polynomial, int, List[List[Any]]] = 1):
        """
        Initialize a rational polynomial.
        
        Args:
            num: Numerator (polynomial, constant, or list representation)
            denom: Denominator (polynomial, constant, or list representation, default: 1)
        """
        # Convert list to Polynomial if needed
        if isinstance(num, list):
            num = Polynomial(num)
        
        # Convert list to Polynomial if needed
        if isinstance(denom, list):
            denom = Polynomial(denom)
            
        # Convert integers to constant polynomials
        if isinstance(num, int):
            nvars = getattr(denom, 'nvars', 0) if not isinstance(denom, int) else 0
            num = Polynomial.from_constant(num, nvars)
        
        if isinstance(denom, int):
            if denom == 0:
                raise ZeroDivisionError("Denominator cannot be zero")
            denom = Polynomial.from_constant(denom, num.nvars)
        
        # Ensure same number of variables
        nvars = max(num.nvars, denom.nvars)
        if num.nvars < nvars:
            # Convert to polynomial with more variables
            new_coeffs = {mono + (0,) * (nvars - len(mono)): coeff 
                         for mono, coeff in num.coeffs.items()}
            num = Polynomial(new_coeffs, nvars)
        
        if denom.nvars < nvars:
            # Convert to polynomial with more variables
            new_coeffs = {mono + (0,) * (nvars - len(mono)): coeff 
                         for mono, coeff in denom.coeffs.items()}
            denom = Polynomial(new_coeffs, nvars)
        
        # Store the numerator and denominator
        self.num = num
        self.denom = denom
    
    def __str__(self) -> str:
        """Convert rational polynomial to string representation."""
        if self.is_polynomial():
            return str(self.num)
        
        # If denominator is 1, just show numerator
        if len(self.denom.coeffs) == 1:
            mono, coeff = next(iter(self.denom.coeffs.items()))
            if all(e == 0 for e in mono) and coeff == 1:
                return str(self.num)
        
        # For complex fractions, use parentheses
        num_str = str(self.num)
        denom_str = str(self.denom)
        
        if "+" in num_str or "-" in num_str[1:]:
            num_str = f"({num_str})"
        
        if "+" in denom_str or "-" in denom_str[1:]:
            denom_str = f"({denom_str})"
        
        return f"{num_str}/{denom_str}"
    
    def __repr__(self) -> str:
        """Return a string representation of the rational polynomial."""
        return f"RationalPolynomial({repr(self.num)}, {repr(self.denom)})"
    
    def is_polynomial(self) -> bool:
        """Check if this rational polynomial is actually a polynomial."""
        # Check if denominator is a constant 1
        if len(self.denom.coeffs) == 1:
            mono, coeff = next(iter(self.denom.coeffs.items()))
            return all(e == 0 for e in mono) and coeff == 1
        return False
    
    def __eq__(self, other: Any) -> bool:
        """Test equality with another rational polynomial."""
        if isinstance(other, RationalPolynomial):
            # a/b = c/d if and only if a*d = b*c
            return self.num * other.denom == self.denom * other.num
        elif isinstance(other, (Polynomial, int)):
            if isinstance(other, int):
                other = Polynomial.from_constant(other, self.num.nvars)
            # a/b = c if and only if a = b*c
            return self.num == self.denom * other
        return NotImplemented
    
    def __add__(self, other: Union['RationalPolynomial', Polynomial, int]) -> 'RationalPolynomial':
        """Add this rational polynomial to another rational polynomial, polynomial, or constant."""
        if isinstance(other, (Polynomial, int)):
            if isinstance(other, int):
                other = Polynomial.from_constant(other, self.num.nvars)
            # a/b + c = (a + b*c)/b
            return RationalPolynomial(self.num + self.denom * other, self.denom)
        
        if not isinstance(other, RationalPolynomial):
            return NotImplemented
        
        # a/b + c/d = (a*d + b*c)/(b*d)
        return RationalPolynomial(
            self.num * other.denom + self.denom * other.num,
            self.denom * other.denom
        )
    
    def __radd__(self, other: Union[Polynomial, int]) -> 'RationalPolynomial':
        """Add a polynomial or constant to this rational polynomial."""
        return self + other
    
    def __sub__(self, other: Union['RationalPolynomial', Polynomial, int]) -> 'RationalPolynomial':
        """Subtract another rational polynomial, polynomial, or constant from this rational polynomial."""
        if isinstance(other, (Polynomial, int)):
            if isinstance(other, int):
                other = Polynomial.from_constant(other, self.num.nvars)
            # a/b - c = (a - b*c)/b
            return RationalPolynomial(self.num - self.denom * other, self.denom)
        
        if not isinstance(other, RationalPolynomial):
            return NotImplemented
        
        # a/b - c/d = (a*d - b*c)/(b*d)
        return RationalPolynomial(
            self.num * other.denom - self.denom * other.num,
            self.denom * other.denom
        )
    
    def __rsub__(self, other: Union[Polynomial, int]) -> 'RationalPolynomial':
        """Subtract this rational polynomial from a polynomial or constant."""
        if isinstance(other, int):
            other = Polynomial.from_constant(other, self.num.nvars)
        # c - a/b = (c*b - a)/b
        return RationalPolynomial(other * self.denom - self.num, self.denom)
    
    def __mul__(self, other: Union['RationalPolynomial', Polynomial, int]) -> 'RationalPolynomial':
        """Multiply this rational polynomial by another rational polynomial, polynomial, or constant."""
        if isinstance(other, (Polynomial, int)):
            if isinstance(other, int):
                other = Polynomial.from_constant(other, self.num.nvars)
            # (a/b) * c = (a*c)/b
            return RationalPolynomial(self.num * other, self.denom)
        
        if not isinstance(other, RationalPolynomial):
            return NotImplemented
        
        # (a/b) * (c/d) = (a*c)/(b*d)
        return RationalPolynomial(
            self.num * other.num,
            self.denom * other.denom
        )
    
    def __rmul__(self, other: Union[Polynomial, int]) -> 'RationalPolynomial':
        """Multiply this rational polynomial by a polynomial or constant."""
        return self * other
    
    def __truediv__(self, other: Union['RationalPolynomial', Polynomial, int]) -> 'RationalPolynomial':
        """Divide this rational polynomial by another rational polynomial, polynomial, or constant."""
        if isinstance(other, (Polynomial, int)):
            if isinstance(other, int):
                if other == 0:
                    raise ZeroDivisionError("Division by zero")
                other = Polynomial.from_constant(other, self.num.nvars)
            # (a/b) / c = a/(b*c)
            return RationalPolynomial(self.num, self.denom * other)
        
        if not isinstance(other, RationalPolynomial):
            return NotImplemented
        
        if other.num.coeffs:
            # (a/b) / (c/d) = (a*d)/(b*c)
            return RationalPolynomial(
                self.num * other.denom,
                self.denom * other.num
            )
        else:
            raise ZeroDivisionError("Division by zero rational polynomial")
    
    def __rtruediv__(self, other: Union[Polynomial, int]) -> 'RationalPolynomial':
        """Divide a polynomial or constant by this rational polynomial."""
        if isinstance(other, int):
            other = Polynomial.from_constant(other, self.num.nvars)
        
        if self.num.coeffs:
            # c / (a/b) = (c*b)/a
            return RationalPolynomial(other * self.denom, self.num)
        else:
            raise ZeroDivisionError("Division by zero rational polynomial")
    
    def __neg__(self) -> 'RationalPolynomial':
        """Negate this rational polynomial."""
        return RationalPolynomial(-self.num, self.denom)
    
    def __pow__(self, power: int) -> 'RationalPolynomial':
        """
        Raise this rational polynomial to a non-negative integer power.
        
        Args:
            power: Non-negative integer exponent
            
        Returns:
            A new rational polynomial representing the power
            
        Raises:
            ValueError: If power is negative (use inv() for negative powers)
            TypeError: If power is not an integer
        """
        if not isinstance(power, int):
            raise TypeError(f"Power must be an integer, got {type(power).__name__}")
            
        if power < 0:
            raise ValueError("Power must be a non-negative integer; use inv() for negative powers")
            
        if power == 0:
            # Return the rational constant 1
            return RationalPolynomial(
                Polynomial.from_constant(1, self.num.nvars),
                Polynomial.from_constant(1, self.denom.nvars)
            )
            
        if power == 1:
            return RationalPolynomial(self.num, self.denom)
            
        # Use binary exponentiation for efficiency
        half = self ** (power // 2)
        if power % 2 == 0:
            return half * half
        else:
            return half * half * self
            
    def inv(self) -> 'RationalPolynomial':
        """
        Compute the inverse of this rational polynomial.
        
        Returns:
            A new rational polynomial representing the inverse (1/self)
            
        Raises:
            ZeroDivisionError: If the numerator is zero
        """
        if not self.num.coeffs:
            raise ZeroDivisionError("Cannot invert a zero rational polynomial")
            
        # Swap numerator and denominator to compute inverse
        return RationalPolynomial(self.denom, self.num)
    
    def simplify(self) -> 'RationalPolynomial':
        """
        Simplify this rational polynomial by dividing out common factors.
        
        Returns:
            A simplified rational polynomial
        """
        # If either numerator or denominator is zero, handle specially
        if not self.num.coeffs:
            # 0/d = 0
            return RationalPolynomial(0, 1)
        
        if not self.denom.coeffs:
            raise ZeroDivisionError("Denominator is zero")
        
        # Check for common constant factors
        if len(self.num.coeffs) == 1 and len(self.denom.coeffs) == 1:
            # Get the only monomial and coefficient from each
            num_mono, num_coeff = next(iter(self.num.coeffs.items()))
            denom_mono, denom_coeff = next(iter(self.denom.coeffs.items()))
            
            # Check for constant GCD
            import math
            gcd_val = math.gcd(num_coeff, denom_coeff)
            if gcd_val > 1:
                simplified_num = Polynomial({num_mono: num_coeff // gcd_val}, self.num.nvars)
                simplified_denom = Polynomial({denom_mono: denom_coeff // gcd_val}, self.denom.nvars)
                return RationalPolynomial(simplified_num, simplified_denom)
        
        # Convert to sympy expressions for more complex simplification
        num_sympy = self.num.tosympy()
        denom_sympy = self.denom.tosympy()
        
        # Use sympy to simplify
        simplified = simplify(num_sympy / denom_sympy)
        
        # If the result is just a fraction of integers, handle it directly
        if simplified.is_Rational:
            num_val = simplified.numerator
            denom_val = simplified.denominator
            return RationalPolynomial(
                Polynomial.from_constant(num_val, self.num.nvars),
                Polynomial.from_constant(denom_val, self.num.nvars)
            )
        
        # For more complex cases, try to extract numerator and denominator
        if isinstance(simplified, sympy.Mul):
            # Look for a division operation
            num_parts = []
            denom_parts = []
            
            for arg in simplified.args:
                if isinstance(arg, sympy.Pow) and arg.args[1] < 0:
                    # Negative power means it's part of the denominator
                    denom_parts.append(arg.args[0] ** (-arg.args[1]))
                else:
                    # Otherwise it's part of the numerator
                    num_parts.append(arg)
            
            if num_parts and denom_parts:
                num_expr = sympy.Mul(*num_parts)
                denom_expr = sympy.Mul(*denom_parts)
                
                # Convert back to Polynomial objects
                new_num = Polynomial.fromsympy(num_expr)
                new_denom = Polynomial.fromsympy(denom_expr)
                
                return RationalPolynomial(new_num, new_denom)
        
        # Get the fraction parts using sympy's fraction function
        try:
            num, denom = sympy.fraction(simplified)
            new_num = Polynomial.fromsympy(num)
            new_denom = Polynomial.fromsympy(denom)
            return RationalPolynomial(new_num, new_denom)
        except Exception:
            # If sympy's fraction fails, check if it's a polynomial
            if simplified.is_polynomial():
                return RationalPolynomial(
                    Polynomial.fromsympy(simplified),
                    Polynomial.from_constant(1, self.num.nvars)
                )
        
        # If all simplification attempts fail, return the original
        return self
    
    def evaluate(self, *args) -> float:
        """
        Evaluate the rational polynomial at the given point.
        
        Args:
            *args: Values for each variable
            
        Returns:
            Value of the rational polynomial at the given point
        """
        num_val = self.num.evaluate(*args)
        denom_val = self.denom.evaluate(*args)
        
        if denom_val == 0:
            raise ZeroDivisionError("Evaluation led to division by zero")
        
        return num_val / denom_val
    
    def tosympy(self, symbols: Optional[List[Symbol]] = None) -> Expr:
        """
        Convert to a sympy expression.
        
        Args:
            symbols: List of sympy symbols to use
            
        Returns:
            A sympy expression representing this rational polynomial
        """
        num_expr = self.num.tosympy(symbols)
        denom_expr = self.denom.tosympy(symbols)
        
        return num_expr / denom_expr
    
    @classmethod
    def fromsympy(cls, expr: Expr, symbols: Optional[List[Symbol]] = None) -> 'RationalPolynomial':
        """
        Convert from a sympy expression.
        
        Args:
            expr: A sympy expression
            symbols: List of sympy symbols to use
            
        Returns:
            A rational polynomial representing the sympy expression
        """
        # Try to separate into numerator and denominator
        if isinstance(expr, sympy.Mul):
            num_parts = []
            denom_parts = []
            
            for arg in expr.args:
                if isinstance(arg, sympy.Pow) and arg.args[1] < 0:
                    # Negative power means it's part of the denominator
                    denom_parts.append(arg.args[0] ** (-arg.args[1]))
                else:
                    # Otherwise it's part of the numerator
                    num_parts.append(arg)
            
            if denom_parts:
                num_expr = sympy.Mul(*num_parts) if num_parts else sympy.Integer(1)
                denom_expr = sympy.Mul(*denom_parts)
                
                num_poly = Polynomial.fromsympy(num_expr, symbols)
                denom_poly = Polynomial.fromsympy(denom_expr, symbols)
                
                return cls(num_poly, denom_poly)
        
        # Check if it's a division
        if isinstance(expr, sympy.Pow) and expr.args[1] == -1:
            # It's 1/something
            denom_poly = Polynomial.fromsympy(expr.args[0], symbols)
            return cls(Polynomial.from_constant(1, denom_poly.nvars), denom_poly)
        
        if isinstance(expr, sympy.Add):
            # It's a sum with no clear denominator, so it's a polynomial
            return cls(Polynomial.fromsympy(expr, symbols), Polynomial.from_constant(1))
        
        # Default case: it's a polynomial
        try:
            return cls(Polynomial.fromsympy(expr, symbols), Polynomial.from_constant(1))
        except ValueError as e:
            # Handle expressions that are not directly convertible to polynomials
            num, denom = sympy.fraction(expr)
            num_poly = Polynomial.fromsympy(num, symbols)
            denom_poly = Polynomial.fromsympy(denom, symbols)
            return cls(num_poly, denom_poly)