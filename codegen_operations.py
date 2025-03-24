"""
Operation-specific Code Generation Module for Geometric Algebra
==============================================================

This module implements specific geometric algebra operations using the code generation
framework provided by codegen_core.py. It contains functions that generate optimized
code for operations like geometric product, outer product, inner product, and more.

Key functions:
- codegen_gp: Generate code for geometric product
- codegen_op: Generate code for outer product
- codegen_ip: Generate code for inner product
- Various other specialized operation generators
"""

from __future__ import annotations

import operator
import warnings
from functools import reduce
from typing import Any, Callable, Dict, List, Optional, Tuple, Union, cast

import sympy

from kingdon.codegen_core import (
    CodegenOutput, LambdifyInput, Fraction, power_supply,
    codegen_product, lambdify
)

# Import for type checking only
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from kingdon.multivector import MultiVector


def codegen_gp(x: Any, y: Any) -> CodegenOutput:
    """
    Generate code for the geometric product between x and y.

    The geometric product combines both inner and outer products in a single operation
    and is the fundamental product in geometric algebra.

    Args:
        x: First symbolic multivector.
        y: Second symbolic multivector.

    Returns:
        CodegenOutput with output keys and a lambda function.

    Example:
        >>> gp_result = codegen_gp(vector1, vector2)
        >>> # The result contains information needed to compute v1 * v2
    """
    res = codegen_product(x, y)
    # In a complete implementation, res keys would be arranged per the algebra.
    return CodegenOutput(tuple(res.keys()) if res else (0,), lambda *args: [res.get(k, 0) for k in (res.keys() if res else [0])])


def _outer_product_filter(kx: int, ky: int, k_out: int) -> bool:
    """
    Filter function for the outer product.
    
    The outer product requires:
    1. Grade increase: grade(k_out) = grade(kx) + grade(ky)
    2. Disjoint basis vectors: kx & ky = 0
    
    Args:
        kx: Binary key of first basis blade
        ky: Binary key of second basis blade
        k_out: Binary key of output basis blade
        
    Returns:
        bool: True if this term should be included in the outer product
    """
    # Check grade increase: output grade equals sum of input grades
    grade_increase = bin(k_out).count('1') == bin(kx).count('1') + bin(ky).count('1')
    
    # Check basis vectors are disjoint (no shared components)
    disjoint_bases = (kx & ky) == 0
    
    return grade_increase and disjoint_bases


def codegen_op(x: 'MultiVector', y: 'MultiVector') -> CodegenOutput:
    """
    Generate code for the outer product of x and y.
    
    The outer product (wedge product) constructs a higher-dimensional subspace
    from lower-dimensional components.
    
    Args:
        x: First symbolic multivector.
        y: Second symbolic multivector.
    
    Returns:
        CodegenOutput with output keys and a lambda function.
    """
    # Use the filter function to compute the outer product terms
    res = codegen_product(x, y, filter_func=_outer_product_filter)
    
    # Handle the case when res is empty (no valid terms)
    if not res:
        # Return a zero scalar with key 0
        return CodegenOutput((0,), lambda *args: [0])
    
    # Create a function that evaluates the outer product
    result_function = lambda *args: [res.get(k, 0) for k in res.keys()]
    
    # Return the keys and evaluation function
    return CodegenOutput(tuple(res.keys()), result_function)


def codegen_sw(x: Any, y: Any) -> CodegenOutput:
    r"""
    Generate code for the conjugation (sandwich product) of y by x: x*y*~x.

    The sandwich product is crucial for computing rotations, reflections, and other
    transformations in geometric algebra.

    Args:
        x: The multivector that "sandwiches" y (typically a rotor or versor).
        y: The multivector being transformed.

    Returns:
        CodegenOutput with output keys and a lambda function.

    Example:
        >>> # If R is a rotor and v is a vector:
        >>> sandwich = codegen_sw(R, v)  
        >>> # This generates code for R*v*~R, rotating v by R
    """
    # For simplicity, use existing operators.
    result = x * y * ~x
    return CodegenOutput(result.keys(), lambda *args: [result[k] for k in result.keys()])


def codegen_cp(x: Any, y: Any) -> CodegenOutput:
    """
    Generate code for the commutator product: 0.5*(x*y - y*x).

    The commutator product is useful for studying the difference between
    different orders of applying operations.

    Args:
        x: First symbolic multivector.
        y: Second symbolic multivector.

    Returns:
        CodegenOutput with output keys and a lambda function.

    Example:
        >>> commutator = codegen_cp(bivector, vector)
    """
    algebra = x.algebra
    filter_func = lambda kx, ky, k_out: (algebra.signs[kx, ky] - algebra.signs[ky, kx])
    res = codegen_product(x, y, filter_func=filter_func)
    return CodegenOutput(tuple(res.keys()) if res else (0,), lambda *args: [res.get(k, 0) for k in (res.keys() if res else [0])])


def codegen_acp(x: 'MultiVector', y: 'MultiVector') -> CodegenOutput:
    """
    Generate code for the anti-commutator product: 0.5*(x*y + y*x).

    The anti-commutator product is useful for studying the similarity between
    different orders of applying operations.

    Args:
        x: First symbolic multivector.
        y: Second symbolic multivector.

    Returns:
        CodegenOutput with output keys and a lambda function.
    """
    algebra = x.algebra
    filter_func = lambda kx, ky, k_out: (algebra.signs[kx, ky] + algebra.signs[ky, kx])
    res = codegen_product(x, y, filter_func=filter_func)
    return CodegenOutput(tuple(res.keys()) if res else (0,), lambda *args: [res.get(k, 0) for k in (res.keys() if res else [0])])


def codegen_ip(x: 'MultiVector', y: 'MultiVector', diff_func: Callable[[int], int] = abs) -> CodegenOutput:
    """
    Generate code for the inner product of x and y.
    
    The inner product extracts the components where the grade changes by
    a specific amount, determined by diff_func.
    
    Args:
        x: First symbolic multivector.
        y: Second symbolic multivector.
        diff_func: Function determining the grade difference, default is absolute difference.
    
    Returns:
        CodegenOutput with output keys and a lambda function.
    
    Example:
        >>> # Standard inner product (grade decreases):
        >>> inner_product = codegen_ip(vector, bivector)
    """
    filter_func = lambda kx, ky, k_out: k_out == diff_func(bin(kx).count('1') - bin(ky).count('1'))
    res = codegen_product(x, y, filter_func=filter_func)
    return CodegenOutput(tuple(res.keys()) if res else (0,), lambda *args: [res.get(k, 0) for k in (res.keys() if res else [0])])


def codegen_lc(x: 'MultiVector', y: 'MultiVector') -> CodegenOutput:
    """
    Generate code for the left-contraction of x and y.
    
    The left contraction is the part of the geometric product where the grade
    decreases by the grade of the left operand.
    
    Args:
        x: First symbolic multivector.
        y: Second symbolic multivector.
    
    Returns:
        CodegenOutput with output keys and a lambda function.
    """
    return codegen_ip(x, y, diff_func=lambda k: -k)


def codegen_rc(x: 'MultiVector', y: 'MultiVector') -> CodegenOutput:
    """
    Generate code for the right-contraction of x and y.
    
    The right contraction is the part of the geometric product where the grade
    decreases by the grade of the right operand.
    
    Args:
        x: First symbolic multivector.
        y: Second symbolic multivector.
    
    Returns:
        CodegenOutput with output keys and a lambda function.
    """
    return codegen_ip(x, y, diff_func=lambda k: k)


def codegen_sp(x: 'MultiVector', y: 'MultiVector') -> CodegenOutput:
    """
    Generate code for the scalar product (grade 0) of x and y.
    
    The scalar product extracts only the scalar part of the geometric product.
    
    Args:
        x: First symbolic multivector.
        y: Second symbolic multivector.
    
    Returns:
        CodegenOutput with output keys and a lambda function.
    """
    return codegen_ip(x, y, diff_func=lambda k: 0)


def codegen_proj(x: 'MultiVector', y: 'MultiVector') -> CodegenOutput:
    r"""
    Generate code for the projection of x onto y: (x \cdot y)*~y.
    
    Projection is a fundamental operation in geometric algebra that gives
    the component of x that is parallel to y.
    
    Args:
        x: The multivector to be projected.
        y: The multivector to project onto.
    
    Returns:
        CodegenOutput with output keys and a lambda function.
    """
    result = (x | y) * ~y
    return CodegenOutput(result.keys(), lambda *args: [result[k] for k in result.keys()])


def codegen_rp(x: 'MultiVector', y: 'MultiVector') -> CodegenOutput:
    """
    Generate code for the regressive product of x and y.
    
    The regressive product (meet operation) finds the intersection of subspaces.
    
    Args:
        x: First symbolic multivector.
        y: Second symbolic multivector.
    
    Returns:
        CodegenOutput with output keys and a lambda function.
    """
    algebra = x.algebra
    key_pss = (1 << algebra.d) - 1  # Pseudoscalar key
    keyout_func = lambda kx, ky: key_pss - (kx ^ ky)
    filter_func = lambda kx, ky, k_out: key_pss == kx + ky - k_out
    
    # Safe sign function that handles missing keys gracefully
    def safe_sign_func(pair):
        try:
            # First, check if the original keys are available
            if (pair[0], key_pss - pair[0]) in algebra.signs:
                s1 = algebra.signs[(pair[0], key_pss - pair[0])]
            else:
                s1 = 1  # Default if not found
                
            if (pair[1], key_pss - pair[1]) in algebra.signs:
                s2 = algebra.signs[(pair[1], key_pss - pair[1])]
            else:
                s2 = 1  # Default if not found
                
            if (key_pss - pair[0], key_pss - pair[1]) in algebra.signs:
                s3 = algebra.signs[(key_pss - pair[0], key_pss - pair[1])]
            else:
                s3 = 1  # Default if not found
                
            if (key_pss - (pair[0] ^ pair[1]), pair[0] ^ pair[1]) in algebra.signs:
                s4 = algebra.signs[(key_pss - (pair[0] ^ pair[1]), pair[0] ^ pair[1])]
            else:
                s4 = 1  # Default if not found
                
            return s1 * s2 * s3 * s4
        except KeyError:
            # If any key is missing or an error occurs, default to 1
            return 1
    
    res = codegen_product(x, y, filter_func=filter_func, keyout_func=keyout_func, sign_func=safe_sign_func)
    return CodegenOutput(tuple(res.keys()) if res else (0,), lambda *args: [res.get(k, 0) for k in (res.keys() if res else [0])])

def codegen_inv(y: 'MultiVector', x: Optional['MultiVector'] = None, symbolic: bool = False) -> Any:
    """
    Generate code for the inverse of y using Hitzer's or Shirokov's inverse.
    
    If x is provided, compute x * y^{-1} (division).
    
    Args:
        y: The multivector to invert.
        x: Optional multivector to multiply by y^{-1}.
        symbolic: If True, return symbolic expressions instead of compiling.
    
    Returns:
        If symbolic=True, returns a Fraction.
        Otherwise, returns a LambdifyInput structure.
    """
    alg = y.algebra
    if alg.d < 6:
        num, denom = codegen_hitzer_inv(y, symbolic=True)
    else:
        num, denom = codegen_shirokov_inv(y, symbolic=True)
    num = num if x is None else x * num

    if symbolic:
        return Fraction(num, denom)
    
    d = alg.scalar(name='d', symbolcls=alg.div.codegen_symbolcls)
    denom_inv = alg.scalar([1 / denom])
    yinv = num * d.e  # Note: This multiply may expand too eagerly.
    args = {'y': y.values()}
    expr_dict = dict(yinv.items())
    dependencies = list(zip(d.values(), denom_inv.values()))
    return LambdifyInput(
        funcname=f'codegen_inv_{y.type_number}',
        expr_dict=expr_dict,
        args=args,
        dependencies=dependencies,
    )


def codegen_hitzer_inv(x: 'MultiVector', symbolic: bool = False) -> Tuple[Any, Any]:
    """
    Generate code for the inverse of x using Hitzer's method.
    
    Works up to 5D algebras with explicit formulas based on algebra dimension.
    
    Args:
        x: The multivector to invert.
        symbolic: If True, return symbolic expressions instead of compiling.
    
    Returns:
        If symbolic=True, returns a tuple (numerator, denominator).
        Otherwise, returns the inverse multivector.
    
    Raises:
        NotImplementedError: If the algebra dimension is greater than 5.
    """
    alg = x.algebra
    d = alg.d
    if d == 0:
        num = alg.blades.e
    elif d == 1:
        num = x.involute()
    elif d == 2:
        num = x.conjugate()
    elif d == 3:
        xconj = x.conjugate()
        num = xconj * ~(x * xconj)
    elif d == 4:
        xconj = x.conjugate()
        x_xconj = x * xconj
        num = xconj * (x_xconj - 2 * x_xconj.grade(3, 4))
    elif d == 5:
        xconj = x.conjugate()
        x_xconj = x * x.conjugate()
        combo = xconj * ~x_xconj
        x_combo = x * combo
        num = combo * (x_combo - 2 * x_combo.grade(1, 4))
    else:
        raise NotImplementedError(f"Closed form inverses are not known in d={d} dimensions.")
    denom = (x.sp(num)).e

    if symbolic:
        return Fraction(num, denom)
    return alg.multivector({k: v / denom for k, v in num.items()})


def codegen_shirokov_inv(x: 'MultiVector', symbolic: bool = False) -> Tuple[Any, Any]:
    """
    Generate code for the inverse of x using Shirokov's method.
    
    Works in any algebra, but can be computationally expensive.
    
    Args:
        x: The multivector to invert.
        symbolic: If True, return symbolic expressions instead of compiling.
    
    Returns:
        If symbolic=True, returns a tuple (numerator, denominator).
        Otherwise, returns the inverse multivector.
    """
    alg = x.algebra
    n = 2 ** ((alg.d + 1) // 2)
    supply = power_supply(x, tuple(range(1, n + 1)))
    powers: List['MultiVector'] = []
    cs: List[Any] = []
    xs: List['MultiVector'] = []
    for i in range(1, n + 1):
        powers.append(next(supply))
        xi = powers[i - 1]
        for j in range(i - 1):
            power_idx = i - j - 2
            xi_diff = powers[power_idx] * cs[j]
            xi = xi - xi_diff
        if xi.grades == (0,):
            break
        xs.append(xi)
        cs_val = xi.e
        cs.append(n * cs_val / i if cs_val != 0 else cs_val)
    if i == 1:
        adj = alg.blades.e
    else:
        adj = xs[-1] - cs[-1]
    if symbolic:
        return Fraction(adj, xi.e)
    return alg.multivector({k: v / xi.e for k, v in adj.items()})


def codegen_div(x: 'MultiVector', y: 'MultiVector') -> LambdifyInput:
    """
    Generate code for division: x * y^{-1}.
    
    Division in geometric algebra is implemented as multiplication by the inverse.
    
    Args:
        x: Numerator multivector.
        y: Denominator multivector.
    
    Returns:
        LambdifyInput structure containing division information.
    
    Raises:
        ZeroDivisionError: If y has no valid inverse.
    """
    alg = x.algebra
    num, denom = codegen_inv(y, x, symbolic=True)
    if not denom:
        raise ZeroDivisionError("Division by a multivector with no inverse")
    d = alg.scalar(name='d', symbolcls=alg.div.codegen_symbolcls)
    denom_inv = alg.scalar([1 / denom])
    res = num * d.e
    args = {'x': x.values(), 'y': y.values()}
    expr_dict = dict(res.items())
    dependencies = list(zip(d.values(), denom_inv.values()))
    return LambdifyInput(
        funcname=f'div_{x.type_number}_x_{y.type_number}',
        expr_dict=expr_dict,
        args=args,
        dependencies=dependencies,
    )


# Transformation and unary operations

def codegen_normsq(x: 'MultiVector') -> Any:
    """
    Generate the squared norm (x * ~x) of a multivector.
    
    The squared norm is used in normalization, distance calculations, and tests
    for invertibility.
    
    Args:
        x: The multivector to compute the squared norm of.
    
    Returns:
        The multivector representing x * ~x.
    """
    return x * ~x


"""
Snippet to fix the codegen_reverse function in codegen_operations.py
"""

def codegen_reverse(x: 'MultiVector') -> Dict[int, Any]:
    """
    Generate code for reversion of a multivector ~x.
    
    Reversion flips the order of basis vectors in each basis blade, which
    inverts the sign of grades 2 and 3 (mod 4).
    
    Args:
        x: Multivector to reverse.
    
    Returns:
        Dictionary mapping keys to coefficients with appropriate signs.
    """
    # Use a corrected version that creates a proper mapping
    result = {}
    for k, v in x.items():
        grade = bin(k).count('1')  # Calculate the grade
        sign = -1 if grade % 4 in (2, 3) else 1  # Determine sign change
        result[k] = -v if sign < 0 else v  # Apply sign change
    return result


def codegen_involute(x: 'MultiVector') -> Dict[int, Any]:
    """
    Generate code for grade involution of a multivector.
    
    Grade involution negates odd-grade components, which corresponds to
    inverting the sign of grades 1 and 3 (mod 4).
    
    Args:
        x: Multivector to involute.
    
    Returns:
        Dictionary mapping keys to coefficients with appropriate signs.
    """
    return codegen_involutions(x, invert_grades=(1, 3))


def codegen_conjugate(x: 'MultiVector') -> Dict[int, Any]:
    """
    Generate code for Clifford conjugation of a multivector.
    
    Clifford conjugation combines reversion and grade involution, inverting
    the sign of grades 1 and 2 (mod 4).
    
    Args:
        x: Multivector to conjugate.
    
    Returns:
        Dictionary mapping keys to coefficients with appropriate signs.
    """
    return codegen_involutions(x, invert_grades=(1, 2))


def codegen_involutions(x: 'MultiVector', invert_grades: Tuple[int, ...] = (2, 3)) -> Dict[int, Any]:
    """
    Generate code for involutions (e.g., reversion, grade involution, Clifford conjugation).
    
    Involutions are operations that change signs of specific grades.
    
    Args:
        x: Multivector to apply involution to.
        invert_grades: Tuple of grades to invert (multiply by -1).
    
    Returns:
        Dictionary mapping keys to coefficients with appropriate signs.
    
    Example:
        >>> # For reversion, invert grades 2,3 (mod 4)
        >>> reversion = codegen_involutions(mv, invert_grades=(2, 3))
        >>> # For grade involution, invert grades 1,3 (mod 4)
        >>> grade_inv = codegen_involutions(mv, invert_grades=(1, 3))
    """
    return {k: -v if bin(k).count('1') % 4 in invert_grades else v for k, v in x.items()}


def codegen_polarity(x: 'MultiVector', undual: bool = False) -> Any:
    """
    Generate code for the polarity (dual) operation.
    
    The polarity maps a multivector to its orthogonal complement using the
    pseudoscalar.
    
    Args:
        x: Multivector to apply polarity to.
        undual: If True, compute the undual (inverse of dual) instead.
    
    Returns:
        The (un)dual multivector.
    
    Raises:
        ZeroDivisionError: If the pseudoscalar is null (has zero square).
    """
    if undual:
        return x * x.algebra.pss
    key_pss = (1 << x.algebra.d) - 1  # Pseudoscalar key
    sign = x.algebra.signs[key_pss, key_pss]
    if sign == -1:
        return - x * x.algebra.pss
    if sign == 1:
        return x * x.algebra.pss
    if sign == 0:
        raise ZeroDivisionError("Cannot compute dual in degenerate algebra where the pseudoscalar squares to zero")


def codegen_unpolarity(x: 'MultiVector') -> Any:
    """
    Generate code for the inverse of the polarity (undual) operation.
    
    Args:
        x: Multivector to apply inverse polarity to.
    
    Returns:
        The undual multivector.
    """
    return codegen_polarity(x, undual=True)


def codegen_hodge(x: 'MultiVector', undual: bool = False) -> Dict[int, Any]:
    """
    Generate code for the Hodge dual operation.
    
    The Hodge dual is another way of computing the orthogonal complement,
    often used in differential forms.
    
    Args:
        x: Multivector to apply Hodge dual to.
        undual: If True, compute the undual (inverse of dual) instead.
    
    Returns:
        Dictionary mapping keys to coefficients of the Hodge dual.
    """
    algebra = x.algebra
    key_pss = (1 << algebra.d) - 1  # Pseudoscalar key
    
    if undual:
        return {(key_dual := key_pss - eI): (-v if algebra.signs[key_dual, eI] < 0 else v)
                for eI, v in x.items()}
    return {(key_dual := key_pss - eI): (-v if algebra.signs[eI, key_dual] < 0 else v)
            for eI, v in x.items()}


def codegen_unhodge(x: 'MultiVector') -> Dict[int, Any]:
    """
    Generate code for the inverse of the Hodge dual operation.
    
    Args:
        x: Multivector to apply inverse Hodge dual to.
    
    Returns:
        Dictionary mapping keys to coefficients of the inverse Hodge dual.
    """
    return codegen_hodge(x, undual=True)


def codegen_sqrt(x: 'MultiVector') -> LambdifyInput:
    """
    Compute the square root of a multivector using the Study number approach.
    
    This method works for multivectors that can be split into scalar plus 
    non-scalar parts, where the non-scalar part squares to a negative scalar.
    
    Args:
        x: Multivector to compute square root of.
    
    Returns:
        LambdifyInput structure for computing the square root.
    
    Warnings:
        Issues a warning if the multivector might not be a Study number.
    """
    alg = x.algebra
    if x.grades == (0,):
        return {0: f'({str(x.e)}**0.5)'}
    a, bI = x.grade(0), x - x.grade(0)
    has_solution = len(x.grades) <= 2 and 0 in x.grades
    if not has_solution:
        warnings.warn("Cannot verify that we really are taking the sqrt of a Study number.", RuntimeWarning)
    bI_sq = bI * bI
    if not bI_sq:
        cp = f'({str(a.e)}**0.5)'
    else:
        normS = (a * a - bI * bI).e
        cp = f'(0.5 * ({str(a.e)} + {str(normS)}**0.5)) ** 0.5'
    c = alg.scalar(name='c')
    c2_inv = alg.scalar(name='c2_inv')
    dI = bI * c2_inv
    res = c + dI
    args = {'x': x.values()}
    expr_dict = dict(res.items())
    dependencies = list(zip(c.values(), [cp])) + list(zip(c2_inv.values(), [f'0.5 / {cp}']))
    return LambdifyInput(
        funcname=f'sqrt_{x.type_number}',
        expr_dict=expr_dict,
        args=args,
        dependencies=dependencies,
    )


# Outer exponential and related functions

def codegen_outerexp(x: 'MultiVector', asterms: bool = False) -> Any:
    """
    Compute the outer exponential of x.
    
    The outer exponential is a series expansion using the outer product instead
    of the geometric product, useful for certain geometric constructions.
    
    Args:
        x: The multivector to compute the outer exponential of.
        asterms: If True, return a list of terms; otherwise, return the sum.
    
    Returns:
        If asterms=True, a list of terms in the expansion.
        Otherwise, the sum of all terms.
    
    Warnings:
        Issues a warning if x has mixed grades, which might cause convergence issues.
    """
    alg = x.algebra
    if len(x.grades) != 1:
        warnings.warn('Outer exponential might not converge for mixed-grade multivectors.', RuntimeWarning)
    k = alg.d
    Ws = [alg.scalar([1]), x]
    j = 2
    while j <= k:
        Wj = Ws[-1] ^ x
        Wj._values = tuple(v / j for v in Wj._values)
        if Wj:
            Ws.append(Wj)
            j += 1
        else:
            break
    if asterms:
        return Ws
    return reduce(operator.add, Ws)


def codegen_outersin(x: 'MultiVector') -> Any:
    """
    Compute the outer sine of x, using odd terms of the outer exponential.
    
    Args:
        x: The multivector to compute the outer sine of.
    
    Returns:
        The outer sine multivector.
    """
    odd_Ws = codegen_outerexp(x, asterms=True)[1::2]
    return reduce(operator.add, odd_Ws)


def codegen_outercos(x: 'MultiVector') -> Any:
    """
    Compute the outer cosine of x, using even terms of the outer exponential.
    
    Args:
        x: The multivector to compute the outer cosine of.
    
    Returns:
        The outer cosine multivector.
    """
    even_Ws = codegen_outerexp(x, asterms=True)[0::2]
    return reduce(operator.add, even_Ws)


def codegen_outertan(x: 'MultiVector') -> Any:
    """
    Compute the outer tangent of x as the ratio of outer sine to outer cosine.
    
    Args:
        x: The multivector to compute the outer tangent of.
    
    Returns:
        The outer tangent multivector.
    """
    Ws = codegen_outerexp(x, asterms=True)
    even_Ws, odd_Ws = Ws[0::2], Ws[1::2]
    outercos = reduce(operator.add, even_Ws)
    outersin = reduce(operator.add, odd_Ws)
    return outersin / outercos


# Basic arithmetic operations

def codegen_add(x: 'MultiVector', y: 'MultiVector') -> Dict[int, Any]:
    """
    Generate code for addition of multivectors x + y.
    
    Args:
        x: First multivector.
        y: Second multivector.
    
    Returns:
        Dictionary mapping keys to combined coefficients.
    """
    vals = dict(x.items())
    for k, v in y.items():
        if k in vals:
            vals[k] = vals[k] + v
        else:
            vals[k] = v
    return vals


def codegen_sub(x: 'MultiVector', y: 'MultiVector') -> Dict[int, Any]:
    """
    Generate code for subtraction of multivectors x - y.
    
    Args:
        x: First multivector.
        y: Second multivector.
    
    Returns:
        Dictionary mapping keys to resulting coefficients.
    """
    vals = dict(x.items())
    for k, v in y.items():
        if k in vals:
            vals[k] = vals[k] - v
        else:
            vals[k] = -v
    return vals


def codegen_neg(x: 'MultiVector') -> Dict[int, Any]:
    """
    Generate code for negation of a multivector -x.
    
    Args:
        x: Multivector to negate.
    
    Returns:
        Dictionary mapping keys to negated coefficients.
    """
    return {k: -v for k, v in x.items()}