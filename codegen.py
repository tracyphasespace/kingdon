"""
Code Generation Module for Geometric Algebra
============================================

This module provides code generation utilities for Geometric Algebra operations.
It creates optimized functions for various products (geometric, outer, inner, etc.)
and other operations on multivectors.

The module allows for efficient computation of Geometric Algebra operations through
code generation, which dynamically creates optimized functions for specific operations.

Example:
    >>> from kingdon.algebra import Algebra
    >>> from kingdon.codegen import do_codegen, codegen_gp
    >>> alg = Algebra(3, 0, 1)  # 3D projective geometric algebra
    >>> x = alg.multivector(name="x", grades=(1,))
    >>> y = alg.multivector(name="y", grades=(1,))
    >>> gp_func = do_codegen(codegen_gp, x, y)
    >>> # gp_func can now be used to compute the geometric product of x and y
"""

# At the top of codegen.py, replace the imports section with:

from __future__ import annotations

import string
import builtins
import keyword
import inspect
import linecache
import operator
import warnings
import math
from collections import Counter, defaultdict, namedtuple
from dataclasses import dataclass
from functools import reduce, cached_property, partial
from itertools import product, combinations, groupby
from typing import (
    Any, Callable, Dict, Generator, Iterable, List, Mapping, 
    NamedTuple, Optional, Set, Tuple, Union, cast, TYPE_CHECKING
)

import numpy as np
from sympy.printing.lambdarepr import LambdaPrinter
from sympy.utilities.iterables import iterable, flatten
import sympy
from sympy import Integer, Symbol

# Defer MultiVector import to function runtime to avoid circular imports
# Use TYPE_CHECKING for type annotations
if TYPE_CHECKING:
    from kingdon.multivector import MultiVector

# Where MultiVector is needed in a function, add a local import:
def codegen_function_example(x: Any) -> Any:
    """Example function showing how to import MultiVector locally."""
    # Local import to avoid circular imports
    from kingdon.multivector import MultiVector
    
    # Function implementation
    # ...
    return x


#=============================================================================
# CORE UTILITIES
#=============================================================================

class mathstr(str):
    """
    Lightweight string subclass that overloads math operators to form symbolic expressions.
    
    This class allows creating mathematical expressions by applying operators to strings,
    handling sign propagation correctly.
    
    Examples:
        >>> a = mathstr('x')
        >>> b = mathstr('y')
        >>> a + b
        'x+y'
        >>> -a
        '-x'
        >>> a * b
        'x*y'
        >>> a - b
        'x-y'
    """
    def __add__(self, other: str) -> mathstr:
        """Add two expressions, handling sign correctly."""
        if other[0] == '-':
            return self.__class__(f'{self}{other}')
        return self.__class__(f'{self}+{other}')

    def __sub__(self, other: str) -> mathstr:
        """Subtract expressions, handling sign correctly."""
        if other[0] == '-':
            return self.__class__(f'{self}+{other[1:]}')
        return self.__class__(f'{self}-{other}')

    def __neg__(self) -> mathstr:
        """Negate an expression, handling double negation correctly."""
        if self[0] == '-':
            return self.__class__(self[1:])
        return self.__class__(f'-{self}')

    def __mul__(self, other: str) -> mathstr:
        """Multiply expressions, handling sign propagation correctly."""
        if other[0] != '-':
            return self.__class__(f'{self}*{other}')
        elif self[0] == '-':
            return self.__class__(f'{self[1:]}*{other[1:]}')
        return self.__class__(f'-{self}*{other[1:]}')


@dataclass
class AdditionChains:
    """
    Computes minimal addition chains for exponentiation with minimal multiplications.
    
    An addition chain for n is a sequence starting with 1 where each element
    is the sum of two previous elements, and the sequence ends with n.
    This allows computing x^n with fewer multiplications than naive exponentiation.
    
    Attributes:
        limit (int): The upper limit for which chains will be computed.
        
    Example:
        >>> chains = AdditionChains(15)
        >>> chains[15]  # Get minimal chain for computing x^15
        (1, 2, 3, 6, 12, 15)
        # To compute x^15, calculate: x^2 = x*x, x^3 = x^2*x, x^6 = x^3*x^3, 
        # x^12 = x^6*x^6, x^15 = x^12*x^3
    """
    limit: int

    @cached_property
    def minimal_chains(self) -> Dict[int, Tuple[int, ...]]:
        """
        Compute minimal addition chains for all integers up to limit.
        
        Returns:
            Dict[int, Tuple[int, ...]]: Dictionary mapping integers to their minimal chains.
        """
        chains: Dict[int, Tuple[int, ...]] = {1: (1,)}
        while any(i not in chains for i in range(1, self.limit + 1)):
            for chain in chains.copy().values():
                right_summand = chain[-1]
                for left_summand in chain:
                    value = left_summand + right_summand
                    if value <= self.limit and value not in chains:
                        chains[value] = (*chain, value)
        return chains

    def __getitem__(self, n: int) -> Tuple[int, ...]:
        """
        Get the minimal addition chain for n.
        
        Args:
            n (int): The target exponent.
            
        Returns:
            Tuple[int, ...]: The minimal addition chain for n.
            
        Raises:
            KeyError: If n is larger than the limit or less than 1.
        """
        return self.minimal_chains[n]

    def __contains__(self, item: Any) -> bool:
        """
        Check if the minimal chain for item has been computed.
        
        Args:
            item (Any): The value to check.
            
        Returns:
            bool: True if the chain exists, False otherwise.
        """
        try:
            _ = self[item]
            return True
        except KeyError:
            return False


class CodegenOutput(NamedTuple):
    """
    Output of a codegen function containing the keys and callable function.
    
    Attributes:
        keys_out (Tuple[int, ...]): Tuple with the output blade binary keys.
        func (Callable[..., List[Any]]): Callable function computing the operation.
    """
    keys_out: Tuple[int, ...]
    func: Callable[..., List[Any]]


class LambdifyInput(NamedTuple):
    """
    Structure for passing input to lambdify.
    
    Attributes:
        funcname (str): The function name.
        args (Dict[str, Any]): A dict of arguments.
        expr_dict (Dict[Any, Any]): Dictionary mapping keys to expressions.
        dependencies (List[Tuple[Any, Any]]): List of dependency (variable, expression) pairs.
    """
    funcname: str
    args: Dict[str, Any]
    expr_dict: Dict[Any, Any]
    dependencies: List[Tuple[Any, Any]]


# Define Fraction as a simple named tuple
Fraction = namedtuple('Fraction', ['numer', 'denom'])
Fraction.__doc__ = "Tuple representing a fraction (numerator and denominator)."


def power_supply(x: Any, exponents: Union[int, Tuple[int, ...]], 
                 operation: Callable[[Any, Any], Any] = operator.mul) -> Generator[Any, None, None]:
    """
    Generate powers of a given multivector using the minimal number of multiplications.

    When called with an int (e.g., power_supply(x, 15)), it yields the sequence of
    multivector powers ending with x^15. When given a tuple of integers, it yields only
    the requested powers.

    Args:
        x: The MultiVector to be raised to powers.
        exponents: Either an int target exponent or a tuple of desired exponents.
        operation: Binary operation to use (default is multiplication).

    Yields:
        Each computed power as a MultiVector.

    Example:
        >>> for p in power_supply(x, (2, 4, 8)):
        ...     print(p)  # Yields x^2, x^4, x^8
    """
    if isinstance(exponents, int):
        target = exponents
        addition_chains = AdditionChains(target)
        exponents = addition_chains[target]
    else:
        addition_chains = AdditionChains(max(exponents))
    
    powers: Dict[int, 'MultiVector'] = {1: x}
    for step in exponents:
        if step not in powers:
            chain = addition_chains[step]
            # Use the second-to-last step in the chain to compute the new power.
            powers[step] = operation(powers[chain[-2]], powers[step - chain[-2]])
        yield powers[step]


def codegen_product(x: Any, y: Any, 
                    filter_func: Optional[Callable[[int, int, int], bool]] = None,
                    sign_func: Optional[Callable[[Tuple[int, int]], int]] = None,
                    keyout_func: Callable[[int, int], int] = operator.xor) -> Dict[int, Any]:
    """
    Helper function for code generation of product-type operations between multivectors.

    This is a general function that supports various geometric algebra products by
    allowing customization of filtering, sign computation, and key combination.

    Args:
        x: First symbolic multivector.
        y: Second symbolic multivector.
        filter_func: Function that returns True if a term should be included, takes (key_x, key_y, key_out).
        sign_func: Function to compute the sign for each term, takes (key_x, key_y). Defaults to using algebra.signs.
        keyout_func: Function to combine binary keys, defaults to XOR (used for geometric product).

    Returns:
        A dictionary mapping output binary keys to symbolic expressions.

    Example:
        >>> # A custom filter for inner product:
        >>> filter_func = lambda kx, ky, k_out: k_out == abs(kx - ky)
        >>> inner_product_terms = codegen_product(x, y, filter_func=filter_func)
    """
    sign_func = sign_func or (lambda pair: x.algebra.signs[pair])
    res: Dict[int, Any] = {}
    
    for (kx, vx), (ky, vy) in product(x.items(), y.items()):
        sign = sign_func((kx, ky))
        if sign:  # Skip if the sign is zero
            key_out = keyout_func(kx, ky)
            if filter_func and not filter_func(kx, ky, key_out):
                continue
            
            termstr = vx * vy if sign > 0 else (-vx * vy)
            if key_out in res:
                res[key_out] += termstr
            else:
                res[key_out] = termstr
                
    return res


def do_codegen(codegen: Callable, *mvs: Any) -> CodegenOutput:
    """
    Generate code for an operation between multivectors.
    
    Args:
        codegen: Function that generates code for the operation
        *mvs: Input multivectors
        
    Returns:
        CodegenOutput with keys and compiled function
    """
    algebra = mvs[0].algebra
    res = codegen(*mvs)
    
    # Handle case where codegen returns a CodegenOutput directly
    if isinstance(res, CodegenOutput):
        return res
    
    # Handle LambdifyInput case
    if isinstance(res, LambdifyInput):
        funcname = res.funcname
        args = res.args
        dependencies = res.dependencies
        res = res.expr_dict
    else:
        # Set default values
        funcname = f'{codegen.__name__}_' + '_x_'.join(f"{mv.type_number}" for mv in mvs)
        args = {arg_name: arg.values() for arg_name, arg in zip(string.ascii_uppercase, mvs)}
        dependencies = None

    # Handle different result types
    if isinstance(res, tuple):  # From codegen_hitzer_inv when not symbolic
        res = res[0]  # Take the multivector part
    
    # Handle scalar, Zero, or similar
    if not hasattr(res, 'keys'):
        if res == 0:
            res = {0: 0}  # Wrap as scalar zero
        else:
            res = {0: res}  # Wrap as scalar
    else:
        # Local import to avoid circular dependencies
        from kingdon.multivector import MultiVector
        
        # Convert MultiVector to dict
        if isinstance(res, MultiVector):
            res = dict(res.items())

    # Sort keys in canonical order if algebra has canon2bin
    if hasattr(algebra, 'canon2bin'):
        # Create a dict with only existing keys matching canon2bin values
        res = {k: res[k] for k in res.keys() if k in algebra.canon2bin.values()}
    
    # Skip CSE for string values
    if not algebra.cse and any(isinstance(v, str) for v in res.values()):
        return func_builder(res, *mvs, funcname=funcname)

    # Extract keys and expressions
    keys, exprs = tuple(res.keys()), list(res.values())
    
    # Use lambdify to create the function
    func = lambdify(args, exprs, funcname=funcname, cse=algebra.cse, dependencies=dependencies)
    
    return CodegenOutput(keys, func)


def func_builder(res_vals: Dict[int, Any], *mvs: Any, funcname: str) -> CodegenOutput:
    """
    Build a Python function for computing products between given multivectors.
    
    Args:
        res_vals: Dictionary mapping keys to expressions.
        *mvs: Input multivectors.
        funcname: Name for the generated function.
    
    Returns:
        CodegenOutput with keys and a callable function.
    """
    args = string.ascii_lowercase[:len(mvs)]
    header = f'def {funcname}({", ".join(args)}):'
    if res_vals:
        body = ''
        for mv, arg in zip(mvs, args):
            body += f'    [{", ".join(str(v) for v in mv.values())}] = {arg}\n'
        return_val = f'    return [{", ".join(str(res_vals.get(k, 0)) for k in res_vals.keys())},]'
    else:
        body = ''
        return_val = '    return list()'
    func_source = f'{header}\n{body}\n{return_val}'
    func_locals = {}
    c = compile(func_source, funcname, 'exec')
    exec(c, {}, func_locals)
    linecache.cache[funcname] = (len(func_source), None, func_source.splitlines(True), funcname)
    func = func_locals[funcname]
    return CodegenOutput(tuple(res_vals.keys()), func)


class KingdonPrinter:
    """
    Custom printer for generating code in the Kingdon library.
    
    This printer is responsible for translating symbolic expressions into
    executable Python code.
    
    Attributes:
        _dummify (bool): Whether to replace symbols with dummy symbols.
        _exprrepr (Callable): Function to convert expressions to strings.
        _argrepr (Callable): Function to convert arguments to strings.
    """
    def __init__(self, printer=None, dummify: bool = False):
        """
        Initialize the KingdonPrinter.
        
        Args:
            printer: Custom printer or function to use for expression representation.
            dummify: Whether to replace symbols with dummy symbols.
        """
        self._dummify = dummify
        from sympy.printing.lambdarepr import LambdaPrinter
        if printer is None:
            printer = LambdaPrinter()
        if inspect.isfunction(printer):
            self._exprrepr = printer
        else:
            if inspect.isclass(printer):
                printer = printer()
            self._exprrepr = printer.doprint
        self._argrepr = LambdaPrinter().doprint

    def doprint(self, funcname: str, args: Any, names: Tuple[str, ...], expr: Any, *, cses=()) -> str:
        """
        Generate Python code for a function with the given expressions.
        
        Args:
            funcname: Name for the generated function.
            args: Arguments for the function.
            names: Names for the arguments.
            expr: Expressions to include in the function.
            cses: Common subexpressions to include.
        
        Returns:
            String containing the generated Python code.
        """
        funcbody: List[str] = []
        if not iterable(args):
            args = [args]
        if cses:
            subvars, subexprs = zip(*cses)
            exprs = [expr] + list(subexprs)
            argstrs, exprs = self._preprocess(args, exprs)
            expr, subexprs = exprs[0], exprs[1:]
            cses = zip(subvars, subexprs)
        else:
            argstrs, expr = self._preprocess(args, expr)
        
        funcargs: List[str] = []
        unpackings: List[str] = []
        for name, argstr in zip(names, argstrs):
            if iterable(argstr):
                funcargs.append(name)
                unpackings.extend(self._print_unpacking(argstr, funcargs[-1]))
            else:
                funcargs.append(argstr)
        
        funcsig = f'def {funcname}({", ".join(funcargs)}):'
        funcbody.extend(self._print_funcargwrapping(funcargs))
        funcbody.extend(unpackings)
        
        for s, e in cses:
            if e is None:
                funcbody.append(f'del {s}')
            else:
                funcbody.append(f'{s} = {self._exprrepr(e)}')
        
        str_expr = _recursive_to_string(self._exprrepr, expr)
        if '\n' in str_expr:
            str_expr = f'({str_expr})'
        funcbody.append(f'return {str_expr}')
        
        funclines = [funcsig] + [f'    {line}' for line in funcbody]
        return '\n'.join(funclines) + '\n'

    @classmethod
    def _is_safe_ident(cls, ident: str) -> bool:
        """
        Check if a string is a valid Python identifier.
        
        Args:
            ident: String to check.
        
        Returns:
            True if the string is a valid Python identifier, False otherwise.
        """
        return isinstance(ident, str) and ident.isidentifier() and not keyword.iskeyword(ident)

    def _preprocess(self, args: Any, expr: Any) -> Tuple[List[str], Any]:
        """
        Preprocess arguments and expressions for code generation.
        
        Args:
            args: Arguments to preprocess.
            expr: Expressions to preprocess.
        
        Returns:
            Tuple of processed argument strings and expressions.
        """
        argstrs = [None] * len(args)
        for i, arg in enumerate(args):
            if iterable(arg):
                s, expr = self._preprocess(arg, expr)
            elif hasattr(arg, 'name'):
                s = arg.name
            elif hasattr(arg, 'is_symbol') and arg.is_symbol:
                s = self._argrepr(arg)
            else:
                s = str(arg)
            argstrs[i] = s
        return argstrs, expr

    def _print_funcargwrapping(self, args: List[str]) -> List[str]:
        """
        Generate code for wrapping function arguments, if needed.
        
        Args:
            args: List of argument names.
        
        Returns:
            List of code lines for argument wrapping.
        """
        return []

    # In the KingdonPrinter class, add or replace the _print_unpacking method:

    def _print_unpacking(self, unpackto: Any, arg: str) -> List[str]:
        """
        Generate code for unpacking arguments into variables.
        
        Args:
            unpackto: Structure to unpack into.
            arg: Argument to unpack.
        
        Returns:
            List of code lines for unpacking.
        """
        if isinstance(unpackto, (list, tuple)):
            # For complex unpackings, use tuple unpacking
            var_names = []
            for i, item in enumerate(unpackto):
                if isinstance(item, str):
                    var_names.append(item)
                else:
                    # Generate a unique variable name
                    var_names.append(f"var_{i}")
            
            # Create the unpacking line
            return [f"{', '.join(var_names)} = {arg}"]
        elif isinstance(unpackto, str):
            # For simple variables, just assign directly
            return [f"{unpackto} = {arg}"]
        else:
            # Default case, use repr to create a valid target name
            return [f"var = {arg}"]

def _recursive_to_string(doprint: Callable, arg: Any) -> str:
    """
    Recursively convert expressions to strings.
    
    Args:
        doprint: Function to convert atomic expressions to strings.
        arg: Expression to convert.
    
    Returns:
        String representation of the expression.
    """
    if isinstance(arg, str):
        return arg
    elif arg is None:
        return 'None'
    elif isinstance(arg, (bool, int, float)):
        return str(arg)
    elif hasattr(arg, 'items') and callable(arg.items):  # Dictionary-like
        items = [f"{_recursive_to_string(doprint, k)}: {_recursive_to_string(doprint, v)}" 
                for k, v in arg.items()]
        return "{" + ", ".join(items) + "}"
    elif hasattr(arg, '__iter__') and not isinstance(arg, (str, sympy.Basic)):  # Iterable but not string/sympy
        if isinstance(arg, list):
            left, right = "[", "]"
        elif isinstance(arg, tuple):
            left, right = "(", ")"
        elif isinstance(arg, set):
            left, right = "{", "}"
        else:
            # Fall back to list representation for unknown iterables
            left, right = "[", "]"
        
        items = [_recursive_to_string(doprint, e) for e in arg]
        if isinstance(arg, tuple) and len(items) == 1:
            # Special case for single-item tuples
            return f"({items[0]},)"
        return f"{left}{', '.join(items)}{right}"
    else:
        # Use the provided printer for other types
        try:
            return doprint(arg)
        except Exception as e:
            # Fall back to str if doprint fails
            return str(arg)


def lambdify(args: Dict[str, Any], exprs: List[Any], funcname: str, 
             dependencies: Optional[List[Tuple[Any, Any]]] = None, 
             printer: Any = LambdaPrinter, dummify: bool = False, cse: bool = False) -> Callable:
    """
    Turn symbolic expressions into a Python function.
    
    This function is inspired by sympy.lambdify but adapted for specialized 
    geometric algebra needs.
    
    Args:
        args: Dictionary mapping parameter names to values.
        exprs: List of expressions to be included in the function.
        funcname: Name for the generated function.
        dependencies: List of (variable, expression) pairs for dependencies.
        printer: Custom printer for code generation.
        dummify: Whether to replace symbols with dummy symbols.
        cse: Whether to use common subexpression elimination.
    
    Returns:
        A callable function that evaluates the expressions.
    """
    if printer is LambdaPrinter:
        printer = LambdaPrinter(
            {'fully_qualified_modules': False, 'inline': True,
             'allow_unknown_functions': True, 'user_functions': {}}
        )
    
    # Convert inputs to sympy expressions where possible
    tosympy = lambda x: x.tosympy() if hasattr(x, 'tosympy') else x
    
    # Handle dictionary inputs specially
    if isinstance(exprs, dict):
        # For dictionary expressions (like from codegen_conjugate),
        # transform into a list of values corresponding to the keys
        keys = sorted(exprs.keys())
        expr_values = [exprs[key] for key in keys]
        # Convert to sympy expressions
        expr_values = [tosympy(e) for e in expr_values]
        # Create return statement that constructs a list
        return_expr = f"[{', '.join(_recursive_to_string(printer.doprint, e) for e in expr_values)}]"
    else:
        # Handle standard list expressions
        args = {name: [tosympy(v) for v in values] for name, values in args.items()}
        exprs = [tosympy(expr) for expr in exprs]
        if dependencies is not None:
            dependencies = [(tosympy(y), tosympy(x)) for y, x in dependencies]
        
        # Process names and arguments
        names = tuple(arg if isinstance(arg, str) else arg.name for arg in args.keys())
        iterable_args = tuple(args.values())
        
        # Use CSE if enabled and expressions are not strings
        lhsides, rhsides = zip(*dependencies) if dependencies else ([], [])
        if cse and not any(isinstance(expr, str) for expr in exprs):
            if not callable(cse):
                from sympy.simplify.cse_main import cse
            if dependencies:
                all_exprs = [*exprs, *rhsides]
                cses, _all_exprs = cse(all_exprs, list=False, order='none', ignore=lhsides)
                _exprs, _rhsides = _all_exprs[:-len(rhsides)], _all_exprs[len(exprs):]
                cses.extend(list(zip(flatten(lhsides), flatten(_rhsides))))
            else:
                cses, _exprs = cse(exprs, list=False)
        else:
            cses, _exprs = list(zip(flatten(lhsides), flatten(rhsides))) if lhsides else [], exprs
        
        if not any(_exprs):
            _exprs = list('0' for _ in _exprs)
        
        # Convert expressions to string representation
        if len(_exprs) == 1:
            return_expr = _recursive_to_string(printer.doprint, _exprs[0])
        else:
            return_expr = f"[{', '.join(_recursive_to_string(printer.doprint, e) for e in _exprs)}]"
    
    # Create the function string
    func_args = ', '.join(str(a) for a in args.keys())
    funcstr = f"def {funcname}({func_args}):\n"
    
    # Add CSE variables if any
    if 'cses' in locals() and cses:
        for s, e in cses:
            e_str = _recursive_to_string(printer.doprint, e)
            funcstr += f"    {s} = {e_str}\n"
    
    # Add the return statement
    if '\n' in return_expr:
        return_expr = f"({return_expr})"
    funcstr += f"    return {return_expr}\n"
    
    # Compile and execute the function
    namespace = {'builtins': builtins, 'range': range}
    funclocals = {}
    filename = f'<{funcname}>'
    c = compile(funcstr, filename, 'exec')
    exec(c, namespace, funclocals)
    linecache.cache[filename] = (len(funcstr), None, funcstr.splitlines(True), filename)
    
    return funclocals[funcname]

#=============================================================================
# OPERATION-SPECIFIC CODE GENERATORS
#=============================================================================

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
    filter_func = lambda kx, ky, k_out: bin(k_out).count('1') == diff_func(bin(kx).count('1') - bin(ky).count('1'))
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

    # Handle zero denominator
    if denom == 0:
        if symbolic:
            return Fraction(num, alg.scalar(0))  # Symbolic representation of undefined inverse
        raise ZeroDivisionError("Multivector has no inverse (denominator is zero)")

    if symbolic:
        return Fraction(num, denom)
    return alg.multivector({k: v / denom for k, v in num.items()}), alg.scalar(1)  # Return num/denom as multivector

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
        sign = -1 if grade % 4 in (2, 3) else 1  # Grades 2 and 3 (mod 4) change sign
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
    
    The Clifford conjugation changes signs for grades 1 and 2 (mod 4).
    
    Args:
        x: Multivector to conjugate
    
    Returns:
        Dictionary mapping keys to conjugated coefficients
    """
    # Import locally to avoid circular import
    from kingdon.multivector import MultiVector
    
    # Create a result dictionary
    result = {}
    
    # Apply conjugation rules to each component
    for k, v in x.items():
        grade = bin(k).count('1')  # Calculate the grade
        sign = -1 if grade % 4 in (1, 2) else 1  # Grades 1 and 2 (mod 4) change sign
        result[k] = sign * v  # Apply sign change
    
    return result

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


# Backward compatibility aliases
# These ensure that existing code can still access the functions
# that were previously spread across multiple modules

# All exports
__all__ = [
    # Core utilities
    'mathstr', 'AdditionChains', 'CodegenOutput', 'LambdifyInput', 'Fraction',
    'power_supply', 'codegen_product', 'do_codegen', 'do_compile', 'func_builder',
    'lambdify', 'KingdonPrinter',
    
    # Product operations
    'codegen_gp', 'codegen_sw', 'codegen_cp', 'codegen_ip', 'codegen_op',
    'codegen_div', 'codegen_rp', 'codegen_acp', 'codegen_proj', 'codegen_sp',
    'codegen_lc', 'codegen_rc',
    
    # Inverse operations
    'codegen_inv', 'codegen_hitzer_inv', 'codegen_shirokov_inv',
    
    # Transformations
    'codegen_normsq', 'codegen_reverse', 'codegen_involute', 'codegen_conjugate',
    'codegen_involutions', 'codegen_sqrt',
    
    # Duality operations
    'codegen_polarity', 'codegen_unpolarity', 'codegen_hodge', 'codegen_unhodge',
    
    # Exponential functions
    'codegen_outerexp', 'codegen_outersin', 'codegen_outercos', 'codegen_outertan',
    
    # Basic arithmetic
    'codegen_add', 'codegen_sub', 'codegen_neg'
]


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