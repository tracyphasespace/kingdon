# Necessary imports (ensure these are at the top of polynomial.py)
from collections import defaultdict
from sympy import Add, Expr, Mul, Symbol, simplify, Integer, sympify, Poly, Pow, S, Number, NumberSymbol, Basic
import warnings
import sympy.polys.polyerrors # For explicit error catching
from typing import Dict, List, Optional, Tuple, Union, Any, Callable, Set, cast # Ensure necessary types are imported

# You might need compare if it's used elsewhere, but it's not used *within* this class
# def compare(m1: Any, m2: Any) -> int: ... (definition omitted for brevity)

# --- Replace the existing Polynomial class with this ---
class Polynomial:
    """
    Represents a multivariate polynomial using SymPy for coefficient handling.

    A polynomial is represented internally as a dictionary mapping monomial tuples
    (representing exponents of variables x0, x1, ...) to SymPy coefficient expressions.
    Example: 3*x0**2*x1 + 2*x0 -> {(2, 1): Integer(3), (1, 0): Integer(2)}

    Attributes:
        coeffs: Dictionary mapping monomial tuples to SymPy coefficients.
        nvars: Number of variables the polynomial is defined over.
    """
    # Type hint for coeffs
    coeffs: Dict[Tuple[int, ...], Expr]

    def __init__(self, coeffs: Optional[Union[Dict[Tuple[int, ...], Any], List[List[Any]]]] = None, nvars: Optional[int] = None):
        """
        Initialize a polynomial.

        Args:
            coeffs: Either a dict mapping monomials to coefficients (int, float, Expr),
                    or a list [[coeff, var1_exp, var2_exp,...], ...] (less common).
                    If None or empty, creates a zero polynomial.
            nvars: Number of variables. If None, inferred from the keys in `coeffs`.
                   Must be provided if `coeffs` is None or empty to set the context.

        Raises:
            TypeError: If coeffs has an invalid format or contains incompatible types.
            ValueError: If inferred nvars conflicts or input is inconsistent.
        """
        processed_coeffs: Dict[Tuple[int, ...], Expr] = {}
        inferred_nvars = 0

        if coeffs is None or (hasattr(coeffs, '__len__') and len(coeffs) == 0):
            # Zero polynomial
            if nvars is None:
                 warnings.warn("Initializing zero Polynomial with nvars=0 as it could not be inferred.", stacklevel=2)
                 self.nvars = 0
            elif isinstance(nvars, int) and nvars >= 0:
                 self.nvars = nvars
            else:
                 raise ValueError("nvars must be a non-negative integer when coeffs is empty.")
            self.coeffs = {} # Explicitly empty dict for zero poly
            return # Added return to ensure no further processing for zero init

        elif isinstance(coeffs, dict):
            if not coeffs: # Empty dict treated as zero polynomial
                 self.coeffs = {}
                 self.nvars = nvars if nvars is not None else 0
                 if self.nvars < 0: raise ValueError("nvars must be non-negative.")
                 return

            max_len = 0
            for mono, coeff in coeffs.items():
                if not isinstance(mono, tuple):
                     raise TypeError(f"Dictionary keys (monomials) must be tuples, got {type(mono).__name__}")
                if not all(isinstance(exp, int) and exp >= 0 for exp in mono):
                     raise ValueError(f"Monomial exponents must be non-negative integers, got {mono}")
                max_len = max(max_len, len(mono))
                try:
                    sympy_coeff = sympify(coeff)
                    if not sympy_coeff.is_zero:
                        processed_coeffs[mono] = sympy_coeff
                except (TypeError, ValueError, sympy.SympifyError) as e:
                     raise TypeError(f"Could not convert coefficient '{coeff}' to SymPy expression: {e}") from e
            inferred_nvars = max_len

            if nvars is None:
                self.nvars = inferred_nvars
            elif isinstance(nvars, int) and nvars >= inferred_nvars:
                self.nvars = nvars
            else:
                 raise ValueError(f"Provided nvars ({nvars}) is less than the length of the longest monomial ({inferred_nvars}).")

            final_coeffs = {}
            padding_zeros = (0,) * (self.nvars - inferred_nvars)
            if padding_zeros: # Only pad if necessary
                for mono, coeff in processed_coeffs.items():
                     final_coeffs[mono + padding_zeros] = coeff
                self.coeffs = final_coeffs
            else:
                 self.coeffs = processed_coeffs


        elif isinstance(coeffs, list):
             if not coeffs: # Empty list treated as zero polynomial
                  self.coeffs = {}
                  self.nvars = nvars if nvars is not None else 0
                  if self.nvars < 0: raise ValueError("nvars must be non-negative.")
                  return

             max_len = 0
             raw_coeffs = {}
             for i, term in enumerate(coeffs):
                 if not isinstance(term, list) or len(term) < 1:
                      raise TypeError(f"Each item in list 'coeffs' must be a non-empty list [coeff, exp0,...], got {term} at index {i}")
                 coeff = term[0]
                 exponents = term[1:]
                 if not all(isinstance(exp, int) and exp >= 0 for exp in exponents):
                      raise ValueError(f"Exponents in list 'coeffs' must be non-negative integers, got {exponents} for term {i}")

                 mono = tuple(exponents)
                 max_len = max(max_len, len(mono))
                 try:
                      sympy_coeff = sympify(coeff)
                      if not sympy_coeff.is_zero:
                           current_coeff = raw_coeffs.get(mono, S.Zero)
                           raw_coeffs[mono] = current_coeff + sympy_coeff
                 except (TypeError, ValueError, sympy.SympifyError) as e:
                      raise TypeError(f"Could not convert coefficient '{coeff}' in list term {i} to SymPy expression: {e}") from e

             inferred_nvars = max_len
             if nvars is None: self.nvars = inferred_nvars
             elif isinstance(nvars, int) and nvars >= inferred_nvars: self.nvars = nvars
             else: raise ValueError(f"Provided nvars ({nvars}) < inferred nvars ({inferred_nvars}) from list.")

             self.coeffs = {}
             vars_to_add = self.nvars - inferred_nvars
             padding = (0,) * vars_to_add
             for mono, coeff in raw_coeffs.items():
                  final_coeff = simplify(coeff) # Simplify accumulated coeff
                  if not final_coeff.is_zero:
                       self.coeffs[mono + padding] = final_coeff

        else:
            raise TypeError(f"Expected dict, list, or None for coeffs, got {type(coeffs).__name__}")

        if not self.coeffs:
             self.coeffs = {}


    @classmethod
    def from_constant(cls, c: Any, nvars: int = 0) -> 'Polynomial':
        """ Create a constant polynomial. """
        if nvars < 0: raise ValueError("nvars must be non-negative.")
        try:
            sympy_c = sympify(c)
            if sympy_c.is_zero: return cls(nvars=nvars) # Zero polynomial
            else:
                mono = (0,) * nvars
                return cls({mono: sympy_c}, nvars=nvars)
        except (TypeError, ValueError, sympy.SympifyError) as e:
             raise TypeError(f"Could not convert constant '{c}' to SymPy expression: {e}") from e

    @classmethod
    def from_variable(cls, var_index: int, nvars: int) -> 'Polynomial':
        """ Create a polynomial representing a single variable x_i. """
        if not isinstance(var_index, int) or not (0 <= var_index < nvars):
            raise ValueError(f"Variable index {var_index} out of range [0, {nvars-1}].")
        if nvars <= 0:
             raise ValueError("nvars must be positive to create a variable.")

        monomial = tuple(1 if i == var_index else 0 for i in range(nvars))
        return cls({monomial: S.One}, nvars=nvars) # Coefficient is SymPy One

    def __str__(self) -> str:
        """Convert polynomial to string representation using SymPy."""
        if not self.coeffs: return "0"
        try:
            return str(self.tosympy())
        except Exception as e:
            warnings.warn(f"Error generating string via SymPy: {e}. Using basic representation.", stacklevel=2)
            items = [f"{coeff if coeff != 1 or all(e == 0 for e in mono) else ''}{'*' if coeff != 1 and any(e > 0 for e in mono) else ''}" +
                     '*'.join(f"x{i}**{exp}" for i, exp in enumerate(mono) if exp > 0)
                     for mono, coeff in self.coeffs.items()]
            return " + ".join(items) if items else "0"


    def __repr__(self) -> str:
        """Return a detailed string representation of the polynomial."""
        sorted_coeffs_repr = dict(sorted(self.coeffs.items()))
        return f"Polynomial({sorted_coeffs_repr}, nvars={self.nvars})"

    def __eq__(self, other: Any) -> bool:
        """Test equality with another polynomial or constant."""
        if isinstance(other, Polynomial):
            nvars = max(self.nvars, other.nvars)
            coeffs1 = self._pad_coeffs(self.coeffs, nvars)
            coeffs2 = other._pad_coeffs(other.coeffs, nvars)
            return coeffs1 == coeffs2
        elif isinstance(other, (int, float, complex, Expr)):
            try:
                other_expr = sympify(other)
                if other_expr.is_zero:
                    return not bool(self.coeffs) # Check if self is zero poly
                if len(self.coeffs) != 1: return False
                mono, coeff = next(iter(self.coeffs.items()))
                if any(e != 0 for e in mono): return False
                return (coeff - other_expr).is_zero
            except (TypeError, ValueError, sympy.SympifyError):
                return False
        return NotImplemented

    def _pad_coeffs(self, coeffs_dict: Dict[Tuple[int,...], Expr], target_nvars: int) -> Dict[Tuple[int,...], Expr]:
        """Helper to pad monomials in a coefficient dictionary to target_nvars."""
        if not coeffs_dict: return {}
        current_nvars = max((len(m) for m in coeffs_dict), default=0)
        vars_to_add = target_nvars - current_nvars
        if vars_to_add < 0:
             raise ValueError(f"Target nvars ({target_nvars}) is less than current nvars ({current_nvars}). Cannot un-pad.")
        if vars_to_add == 0: return coeffs_dict # No padding needed
        padding = (0,) * vars_to_add
        padded = {mono + padding: coeff for mono, coeff in coeffs_dict.items()}
        return padded

    def _align_other(self, other: Union['Polynomial', int, float, complex, Expr]) -> Tuple['Polynomial', 'Polynomial', int]:
        """Align self and other to have the same nvars."""
        if isinstance(other, Polynomial):
            nvars = max(self.nvars, other.nvars)
            self_coeffs = self._pad_coeffs(self.coeffs, nvars)
            other_coeffs = other._pad_coeffs(other.coeffs, nvars)
            aligned_self = Polynomial(self_coeffs, nvars) if nvars > self.nvars else self
            aligned_other = Polynomial(other_coeffs, nvars) if nvars > other.nvars else other
            return aligned_self, aligned_other, nvars
        elif isinstance(other, (int, float, complex, Expr)):
            other_poly = Polynomial.from_constant(other, self.nvars)
            return self, other_poly, self.nvars
        else:
            raise TypeError(f"Unsupported operand type(s) for polynomial operation: Polynomial and {type(other).__name__}")


    def __add__(self, other: Union['Polynomial', int, float, complex, Expr]) -> 'Polynomial':
        """Add this polynomial to another polynomial or constant."""
        try:
            aligned_self, aligned_other, nvars = self._align_other(other)
        except TypeError:
            return NotImplemented
        result_coeffs: Dict[Tuple[int,...], Expr] = defaultdict(lambda: S.Zero)
        for mono, coeff in aligned_self.coeffs.items():
            result_coeffs[mono] += coeff
        for mono, coeff in aligned_other.coeffs.items():
            result_coeffs[mono] += coeff
        final_coeffs = {mono: simp_coeff for mono, coeff in result_coeffs.items() if not (simp_coeff := simplify(coeff)).is_zero}
        return Polynomial(final_coeffs, nvars)

    def __radd__(self, other: Union[int, float, complex, Expr]) -> 'Polynomial':
        """Add a constant/expression to this polynomial."""
        return self + other

    def __sub__(self, other: Union['Polynomial', int, float, complex, Expr]) -> 'Polynomial':
        """Subtract another polynomial or constant."""
        return self + (-other)

    def __rsub__(self, other: Union[int, float, complex, Expr]) -> 'Polynomial':
        """Subtract this polynomial from a constant/expression."""
        try:
            other_poly = Polynomial.from_constant(other, self.nvars)
            return other_poly - self
        except TypeError:
             return NotImplemented


    def __mul__(self, other: Union['Polynomial', int, float, complex, Expr]) -> 'Polynomial':
        """Multiply by another polynomial or constant/expression."""
        try:
             aligned_self, aligned_other, nvars = self._align_other(other)
        except TypeError:
             return NotImplemented
        result_coeffs: Dict[Tuple[int,...], Expr] = defaultdict(lambda: S.Zero)
        for mono1, coeff1 in aligned_self.coeffs.items():
            for mono2, coeff2 in aligned_other.coeffs.items():
                product_coeff = coeff1 * coeff2
                if not product_coeff.is_zero:
                    product_mono = tuple(e1 + e2 for e1, e2 in zip(mono1, mono2))
                    result_coeffs[product_mono] += product_coeff
        final_coeffs = {mono: simp_coeff for mono, coeff in result_coeffs.items() if not (simp_coeff := simplify(coeff)).is_zero}
        return Polynomial(final_coeffs, nvars)


    def __rmul__(self, other: Union[int, float, complex, Expr]) -> 'Polynomial':
        """Multiply this polynomial by a constant/expression."""
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
        if power == 0: return Polynomial.from_constant(1, self.nvars)
        if power == 1: return Polynomial(self.coeffs.copy(), self.nvars) # Return a copy
        result = Polynomial.from_constant(1, self.nvars)
        base = self
        p = power
        while p > 0:
             if p % 2 == 1: result = result * base
             base = base * base
             p //= 2
        return result

    def evaluate(self, *args: Any) -> Union[int, float, complex, Expr]:
        """
        Evaluate the polynomial at the given point (numeric or symbolic).
        """
        if len(args) != self.nvars:
            raise ValueError(f"Expected {self.nvars} values for evaluation, got {len(args)}")
        try:
             sympy_args = [sympify(arg) for arg in args]
        except (TypeError, ValueError, sympy.SympifyError) as e:
             raise TypeError(f"Could not convert evaluation arguments {args} to SymPy expressions: {e}") from e
        result: Expr = S.Zero
        for mono, coeff in self.coeffs.items():
            term: Expr = coeff
            for exp, val in zip(mono, sympy_args):
                if exp > 0:
                     term *= (val ** exp)
            result += term
        return simplify(result)

    def derive(self, var_index: int) -> 'Polynomial':
        """ Compute the partial derivative with respect to variable x_i. """
        if not isinstance(var_index, int) or not (0 <= var_index < self.nvars):
            raise ValueError(f"Variable index {var_index} out of range [0, {self.nvars-1}].")
        result_coeffs: Dict[Tuple[int,...], Expr] = defaultdict(lambda: S.Zero)
        for mono, coeff in self.coeffs.items():
            original_exponent = mono[var_index]
            if original_exponent > 0:
                new_coeff = coeff * original_exponent
                if not new_coeff.is_zero:
                     new_mono_list = list(mono)
                     new_mono_list[var_index] -= 1
                     result_coeffs[tuple(new_mono_list)] += new_coeff
        final_coeffs = {mono: simp_coeff for mono, coeff in result_coeffs.items() if not (simp_coeff := simplify(coeff)).is_zero}
        return Polynomial(final_coeffs, self.nvars)

    def tosympy(self, symbols: Optional[List[Symbol]] = None) -> Expr:
        """ Convert to a SymPy expression. """
        if not self.coeffs: return S.Zero
        if symbols is None:
            symbols = [Symbol(f"x{i}") for i in range(self.nvars)]
        elif len(symbols) != self.nvars:
            raise ValueError(f"Expected {self.nvars} symbols, got {len(symbols)}")
        terms = []
        for mono, coeff in self.coeffs.items():
            term: Expr = coeff
            for i, exp in enumerate(mono):
                if exp > 0:
                    term *= symbols[i] ** exp
            terms.append(term)
        return Add(*terms, evaluate=False).doit()

    @classmethod
    def fromsympy(cls, expr: Expr, symbols: Optional[List[Symbol]] = None) -> 'Polynomial':
        """ Convert from a SymPy expression. """
        try:
             # Ensure expr is a SymPy expression - DO THIS ONCE
             expr = sympify(expr)
        except (TypeError, ValueError, sympy.SympifyError) as e:
             raise TypeError(f"Input 'expr' is not a valid SymPy expression or convertible: {e}") from e

        # Determine and validate symbols - DO THIS ONCE
        if symbols is None:
            symbols = sorted(list(expr.free_symbols), key=str)
        else:
            # Validate provided symbols
            if not all(isinstance(s, Symbol) for s in symbols):
                raise TypeError("Provided 'symbols' must be a list of SymPy Symbols.")

        nvars = len(symbols)

        # --- Debug Print (can be removed later) ---
        # print(f"DEBUG: Checking condition. expr={repr(expr)}, type(expr)={type(expr)}, expr.is_constant={expr.is_constant()}, symbols={symbols}")

        # --- FIX: Handle ANY constant expression when symbols is empty ---
        # Use is_constant() which covers numbers and symbols like S.Pi
        # Use is_number if only numeric constants should be handled this way
        if expr.is_constant() and not symbols:
             # print("DEBUG: Caught constant expr and empty symbols. Returning constant polynomial.") # Optional debug
             # Check if it's zero, otherwise return constant poly {(): expr}
             if expr.is_zero:
                 return cls({}, nvars=0) # Zero polynomial
             else:
                 # Create polynomial representing the constant `expr` with nvars=0
                 # The key () represents the constant term (all variables to power 0)
                 return cls({(): expr}, nvars=0)

        # --- Debug Print (can be removed later) ---
        # print("DEBUG: Did NOT catch constant/empty case. Proceeding to Poly().")

        # --- Call Poly() - DO THIS ONCE ---
        try:
            # Convert the expression to a SymPy Poly object w.r.t the symbols
            poly_obj = Poly(expr, *symbols)
        except sympy.PolynomialError as e:
            # Handle cases where expr is not polynomial in the given symbols
            raise ValueError(f"Input expression '{expr}' is not a polynomial in symbols {symbols}. Error: {e}") from e
        except sympy.polys.polyerrors.GeneratorsNeeded as e_gen:
            # Catch GeneratorsNeeded specifically
            # This should now only happen if expr is non-constant AND symbols is empty,
            # or other more complex SymPy internal issues.
            raise ValueError(f"Error converting expression '{expr}' to Poly object (symbols={symbols}): {e_gen}. Ensure generators are provided if needed.") from e_gen
        except Exception as e: # Catch other potential errors during Poly creation
            raise ValueError(f"Error converting expression '{expr}' to Poly object (symbols={symbols}): {e}") from e

        # --- Process results - DO THIS ONCE ---
        coeffs_dict = {}
        for mono_exponents, coeff in poly_obj.terms():
             simp_coeff = simplify(coeff)
             if not simp_coeff.is_zero:
                  # Ensure monomial tuple has length nvars (Poly should handle this)
                  if len(mono_exponents) != nvars:
                       raise RuntimeError(f"Internal Error: SymPy Poly monomial length {len(mono_exponents)} != nvars {nvars}.")
                  coeffs_dict[mono_exponents] = simp_coeff

        return cls(coeffs_dict, nvars=nvars) if coeffs_dict else cls({}, nvars=nvars)


# --- RationalPolynomial Class Definition ---
# (Rest of the file including RationalPolynomial class remains unchanged)
# ... (Previous RationalPolynomial code goes here) ...

class RationalPolynomial:
    """ Represents a rational polynomial (fraction of polynomials). """
    num: Polynomial
    denom: Polynomial

    def __init__(self, num: Union[Polynomial, int, float, complex, Expr, List[List[Any]]],
                 denom: Union[Polynomial, int, float, complex, Expr, List[List[Any]]] = 1):
        """
        Initialize a rational polynomial P/Q.

        Args:
            num: Numerator (Polynomial, constant, SymPy Expr, or list format).
            denom: Denominator (default: 1). Cannot be zero.

        Raises:
            TypeError: If num or denom are invalid types.
            ZeroDivisionError: If the denominator simplifies to zero.
            ValueError: If nvars inferred from num and denom conflict.
        """
        # --- Convert Numerator and Denominator to Polynomial objects ---
        try:
            if isinstance(num, Polynomial): num_poly = num
            # Use Polynomial.from_constant for simple numeric types for efficiency/clarity
            elif isinstance(num, (int, float, complex)): num_poly = Polynomial.from_constant(num)
            else: num_poly = Polynomial.fromsympy(sympify(num)) # Use fromsympy for robustness on expressions/lists
        except (TypeError, ValueError, sympy.SympifyError) as e:
            raise TypeError(f"Invalid numerator input '{num}': {e}") from e

        try:
            if isinstance(denom, Polynomial): denom_poly = denom
            elif isinstance(denom, (int, float, complex)): denom_poly = Polynomial.from_constant(denom)
            else: denom_poly = Polynomial.fromsympy(sympify(denom))
        except (TypeError, ValueError, sympy.SympifyError) as e:
            raise TypeError(f"Invalid denominator input '{denom}': {e}") from e

        # --- Check for Zero Denominator ---
        if not denom_poly.coeffs: # Check if denominator polynomial is zero
            # Check specifically if it was constant 0
            if isinstance(denom, (int, float, complex)) and denom == 0:
                 raise ZeroDivisionError("Denominator cannot be zero.")
            # Or if a more complex expression simplified to zero
            warnings.warn(f"Denominator expression '{denom}' simplified to zero polynomial.", stacklevel=2)
            raise ZeroDivisionError("Denominator polynomial simplifies to zero.")

        # --- Align nvars ---
        nvars = max(num_poly.nvars, denom_poly.nvars)
        # Need to handle alignment correctly, especially if one is constant (nvars=0 initially)
        self.num = Polynomial(num_poly._pad_coeffs(num_poly.coeffs, nvars), nvars) if num_poly.nvars < nvars or not num_poly.coeffs else num_poly
        self.denom = Polynomial(denom_poly._pad_coeffs(denom_poly.coeffs, nvars), nvars) if denom_poly.nvars < nvars or not denom_poly.coeffs else denom_poly


    def __str__(self) -> str:
        """Convert rational polynomial to string representation."""
        # Use SymPy for potentially simplified string representation
        try: return str(self.tosympy())
        except Exception: pass # Fallback if tosympy fails

        # Fallback formatting
        if self.is_polynomial(): return str(self.num)
        num_str = str(self.num); denom_str = str(self.denom)
        # Add parentheses if needed (check if result of tosympy is Add)
        num_sympy = self.num.tosympy()
        den_sympy = self.denom.tosympy()
        num_has_terms = isinstance(num_sympy, Add)
        den_has_terms = isinstance(den_sympy, Add)
        if num_has_terms: num_str = f"({num_str})"
        if den_has_terms: denom_str = f"({denom_str})"
        # Avoid unnecessary /1
        if denom_str == "1": return num_str
        return f"{num_str}/{denom_str}"

    def __repr__(self) -> str:
        """Return a detailed string representation."""
        return f"RationalPolynomial(num={repr(self.num)}, denom={repr(self.denom)})"

    def is_polynomial(self) -> bool:
        """Check if denominator is equivalent to Polynomial 1."""
        one_poly = Polynomial.from_constant(1, self.denom.nvars)
        return self.denom == one_poly

    def __eq__(self, other: Any) -> bool:
        """Test equality: self == other"""
        if isinstance(other, RationalPolynomial):
            # a/b == c/d <=> a*d == b*c
            # Align nvars for comparison
            nvars = max(self.num.nvars, other.num.nvars)
            num1 = Polynomial(self.num._pad_coeffs(self.num.coeffs, nvars), nvars)
            den1 = Polynomial(self.denom._pad_coeffs(self.denom.coeffs, nvars), nvars)
            num2 = Polynomial(other.num._pad_coeffs(other.num.coeffs, nvars), nvars)
            den2 = Polynomial(other.denom._pad_coeffs(other.denom.coeffs, nvars), nvars)
            try:
                 # Perform cross-multiplication and check equality
                 return (num1 * den2) == (den1 * num2)
            except Exception: return False # Error during operation means not equal
        elif isinstance(other, (Polynomial, int, float, complex, Expr)):
            # Convert other to RationalPolynomial for comparison
            try: other_rp = RationalPolynomial(other, 1)
            except (TypeError, ValueError): return False # Cannot convert other
            return self == other_rp # Reuse RationalPolynomial equality
        return NotImplemented

    def _align_other_rp(self, other: Union['RationalPolynomial', Polynomial, int, float, complex, Expr]) -> Tuple['RationalPolynomial', 'RationalPolynomial']:
        """Convert other to RationalPolynomial and align nvars with self."""
        if isinstance(other, RationalPolynomial): other_rp = other
        elif isinstance(other, (Polynomial, int, float, complex, Expr)):
            try: other_rp = RationalPolynomial(other, 1)
            except (TypeError, ValueError) as e: raise TypeError(f"Cannot convert operand {type(other).__name__} to RationalPolynomial: {e}") from e
        else: raise TypeError(f"Unsupported operand type for RationalPolynomial operation: {type(other).__name__}")

        # Align nvars
        nvars = max(self.num.nvars, other_rp.num.nvars)
        aligned_self_num = Polynomial(self.num._pad_coeffs(self.num.coeffs, nvars), nvars)
        aligned_self_den = Polynomial(self.denom._pad_coeffs(self.denom.coeffs, nvars), nvars)
        aligned_other_num = Polynomial(other_rp.num._pad_coeffs(other_rp.num.coeffs, nvars), nvars)
        aligned_other_den = Polynomial(other_rp.denom._pad_coeffs(other_rp.denom.coeffs, nvars), nvars)

        # Return aligned operands as new objects
        aligned_self = RationalPolynomial(aligned_self_num, aligned_self_den)
        aligned_other = RationalPolynomial(aligned_other_num, aligned_other_den)
        return aligned_self, aligned_other

    def __add__(self, other: Union['RationalPolynomial', Polynomial, int, float, complex, Expr]) -> 'RationalPolynomial':
        """Add: self + other"""
        try: aligned_self, aligned_other = self._align_other_rp(other)
        except TypeError: return NotImplemented

        # a/b + c/d = (a*d + b*c) / (b*d)
        new_num = aligned_self.num * aligned_other.denom + aligned_self.denom * aligned_other.num
        new_den = aligned_self.denom * aligned_other.denom
        # Optional: simplify here or let user call simplify()
        return RationalPolynomial(new_num, new_den).simplify() # Simplify result

    def __radd__(self, other: Union[Polynomial, int, float, complex, Expr]) -> 'RationalPolynomial':
        """Add: other + self"""
        # Addition is commutative
        return self + other

    def __sub__(self, other: Union['RationalPolynomial', Polynomial, int, float, complex, Expr]) -> 'RationalPolynomial':
        """Subtract: self - other"""
        return (self + (-other)).simplify() # Reuse add and neg, simplify result

    def __rsub__(self, other: Union[Polynomial, int, float, complex, Expr]) -> 'RationalPolynomial':
        """Subtract: other - self"""
        try: other_rp, aligned_self = self._align_other_rp(other) # Align self to other's potential nvars
        except TypeError: return NotImplemented
        return (other_rp - aligned_self).simplify() # Reuse sub, simplify result

    def __mul__(self, other: Union['RationalPolynomial', Polynomial, int, float, complex, Expr]) -> 'RationalPolynomial':
        """Multiply: self * other"""
        try: aligned_self, aligned_other = self._align_other_rp(other)
        except TypeError: return NotImplemented

        # (a/b) * (c/d) = (a*c) / (b*d)
        new_num = aligned_self.num * aligned_other.num
        new_den = aligned_self.denom * aligned_other.denom
        return RationalPolynomial(new_num, new_den).simplify() # Simplify result

    def __rmul__(self, other: Union[Polynomial, int, float, complex, Expr]) -> 'RationalPolynomial':
        """Multiply: other * self"""
        # Multiplication is commutative
        return self * other

    def __truediv__(self, other: Union['RationalPolynomial', Polynomial, int, float, complex, Expr]) -> 'RationalPolynomial':
        """Divide: self / other"""
        try: aligned_self, aligned_other = self._align_other_rp(other)
        except TypeError: return NotImplemented

        # Check for division by zero rational polynomial
        if not aligned_other.num.coeffs:
            raise ZeroDivisionError("Division by zero RationalPolynomial.")

        # (a/b) / (c/d) = (a*d) / (b*c)
        new_num = aligned_self.num * aligned_other.denom
        new_den = aligned_self.denom * aligned_other.num

        # Check if the new denominator is zero (can happen with symbolic cancellation)
        if not new_den.coeffs:
             # Check if original denominator was non-zero constant that got cancelled
             if aligned_self.denom.coeffs and not aligned_other.num.coeffs:
                  # Dividing non-zero by zero
                   raise ZeroDivisionError("Division resulted in division by zero.")
             # Otherwise, likely 0/0 or similar undefined scenario
             warnings.warn("Division resulted in an undefined form (potentially 0/0).", stacklevel=2)
             # Returning 0 or raising might be appropriate depending on desired semantics
             # Let's return 0 for now, assuming non-problematic cancellation led to 0/non-zero
             return RationalPolynomial(0, 1) # Zero rational polynomial


        return RationalPolynomial(new_num, new_den).simplify() # Simplify result

    def __rtruediv__(self, other: Union[Polynomial, int, float, complex, Expr]) -> 'RationalPolynomial':
        """Divide: other / self"""
        try: other_rp, aligned_self = self._align_other_rp(other) # Align self
        except TypeError: return NotImplemented

        # other / self = other * self.inv()
        return (other_rp * aligned_self.inv()).simplify() # Reuse mul and inv, simplify

    def __neg__(self) -> 'RationalPolynomial':
        """Negate: -self"""
        # Create new object to avoid modifying self
        return RationalPolynomial(-self.num, self.denom)

    def __pow__(self, power: int) -> 'RationalPolynomial':
        """Power: self ** power"""
        if not isinstance(power, int):
            raise TypeError(f"Power must be an integer, got {type(power).__name__}")

        if power == 0: return RationalPolynomial(1, 1) # Rational 1
        if power == 1: return RationalPolynomial(self.num, self.denom) # Return copy

        if power < 0:
            # Handle negative power: (a/b)^-n = (b/a)^n
            inv_self = self.inv() # Let inv handle zero numerator error
            return inv_self ** (-power) # Simplify handled by recursive call

        # Use Polynomial power method, simplify at the end
        new_num = self.num ** power
        new_den = self.denom ** power
        return RationalPolynomial(new_num, new_den).simplify()


    def inv(self) -> 'RationalPolynomial':
        """Compute the inverse (1/self = denom/num)."""
        if not self.num.coeffs: # Check if numerator polynomial is zero
            raise ZeroDivisionError("Cannot invert a zero RationalPolynomial (zero numerator).")
        # Swap numerator and denominator
        # Return a new instance, do not simplify here as inv might be intermediate step
        return RationalPolynomial(self.denom, self.num)

    def simplify(self) -> 'RationalPolynomial':
        """ Simplify using SymPy's cancel function. """
        try:
            # Convert to a single SymPy expression, cancel, then convert back
            nvars = self.num.nvars
            symbols = [Symbol(f"x{i}") for i in range(nvars)]
            sympy_expr = self.tosympy(symbols)

            # Use sympy.cancel to simplify the rational function
            simplified_expr = sympy.cancel(sympy_expr)

            # Check if simplification resulted in non-rational (e.g., sqrt)
            if not simplified_expr.is_rational_function(*symbols):
                 warnings.warn(f"Simplification resulted in non-rational expression: {simplified_expr}. Returning original.", stacklevel=2)
                 return RationalPolynomial(self.num, self.denom) # Return original

            # Convert back to RationalPolynomial using its numerator/denominator
            num_simp, den_simp = simplified_expr.as_numer_denom()
            # Convert SymPy num/den back to Polynomial objects
            num_poly = Polynomial.fromsympy(num_simp, symbols)
            den_poly = Polynomial.fromsympy(den_simp, symbols)
            return RationalPolynomial(num_poly, den_poly) # Already simplified

        except Exception as e:
            warnings.warn(f"SymPy simplification failed: {e}. Returning original.", stacklevel=2)
            return RationalPolynomial(self.num, self.denom) # Return original copy

    def evaluate(self, *args: Any) -> Union[int, float, complex, Expr]:
        """ Evaluate at the given point (numeric or symbolic). """
        # Ensure correct number of arguments matches nvars
        nvars = self.num.nvars # Num and denom guaranteed to have same nvars
        if len(args) != nvars:
            raise ValueError(f"Expected {nvars} values for evaluation, got {len(args)}")

        try:
             num_val = self.num.evaluate(*args)
             denom_val = self.denom.evaluate(*args)
        except Exception as e:
             raise ValueError(f"Error evaluating numerator/denominator polynomials: {e}") from e

        # Check for division by zero (numeric or symbolic)
        is_zero_denom = False
        if isinstance(denom_val, (int, float, complex, Number, NumberSymbol)):
             # Use tighter tolerance for zero check after evaluation
             is_zero_denom = abs(complex(denom_val)) < 1e-15
        elif hasattr(denom_val, 'is_zero'): # SymPy expression
             try: is_zero_denom = denom_val.is_zero
             except TypeError: pass # is_zero might not be computable

        if is_zero_denom:
            raise ZeroDivisionError(f"Evaluation led to division by zero (denominator evaluated to {denom_val}).")

        # Perform division (SymPy handles symbolic/numeric division)
        try:
             # Sympify ensure both are sympy objects before division
             result_expr = sympify(num_val) / sympify(denom_val)
             # Try simplifying the final result
             return simplify(result_expr)
        except Exception as e:
             raise ValueError(f"Error during final division/simplification in evaluation: {e}") from e


    def tosympy(self, symbols: Optional[List[Symbol]] = None) -> Expr:
        """ Convert to a SymPy expression P/Q. """
        nvars = self.num.nvars
        if symbols is None: symbols = [Symbol(f"x{i}") for i in range(nvars)]
        elif len(symbols) != nvars: raise ValueError(f"Expected {nvars} symbols, got {len(symbols)}")

        num_expr = self.num.tosympy(symbols)
        denom_expr = self.denom.tosympy(symbols)

        # Avoid division if denominator is 1
        if denom_expr == S.One:
            return num_expr
        # Return symbolic fraction, do not simplify here
        return num_expr / denom_expr

    @classmethod
    def fromsympy(cls, expr: Expr, symbols: Optional[List[Symbol]] = None) -> 'RationalPolynomial':
        """ Convert from a SymPy expression P or P/Q. """
        # NOTE: This method on RationalPolynomial might need adjustment if called directly.
        # Currently, it just creates a Polynomial, which might not be the intended behavior
        # if expr is truly rational (P/Q).
        # For now, assuming it behaves like Polynomial.fromsympy for compatibility,
        # but ideally it should handle rational expressions properly.
        # It reuses Polynomial.fromsympy logic by calling it below.

        # Let's delegate to Polynomial.fromsympy for now and wrap the result
        # This assumes the caller wants a RationalPolynomial P/1 if expr is polynomial
        try:
             poly_num = Polynomial.fromsympy(expr, symbols)
             return cls(poly_num, 1) # Create RationalPolynomial P/1
        except ValueError as e_poly:
             # If Polynomial.fromsympy fails (e.g., expr not polynomial),
             # maybe try parsing as num/den directly?
             try:
                  expr_s = sympify(expr)
                  if symbols is None: symbols = sorted(list(expr_s.free_symbols), key=str)
                  num_ex, den_ex = expr_s.as_numer_denom()
                  num_p = Polynomial.fromsympy(num_ex, symbols)
                  den_p = Polynomial.fromsympy(den_ex, symbols)
                  return cls(num_p, den_p)
             except Exception as e_rational:
                  raise ValueError(f"Input expression '{expr}' could not be converted to Polynomial or parsed as RationalPolynomial: {e_poly}; {e_rational}") from e_rational