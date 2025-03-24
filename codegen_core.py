"""
Core Code Generation Module for Geometric Algebra
================================================

This module provides the foundation and utilities for code generation in Geometric Algebra.
It includes base classes, helper functions, and core infrastructure for generating
optimized code for various operations.

Key components:
- CodegenOutput: Container for output keys and function
- LambdifyInput: Structure for passing input to lambdify
- AdditionChains: Utility for computing powers with minimal multiplications
- mathstr: String subclass for symbolic expression building
- Helper functions for code generation and compilation
"""

from __future__ import annotations

import string
import builtins
import keyword
import inspect
import linecache
import operator
import warnings
from collections import defaultdict, namedtuple
from dataclasses import dataclass
from functools import reduce, cached_property
from itertools import product, combinations, groupby
from typing import (
    Any, Callable, Dict, Generator, Iterable, List, Mapping, 
    NamedTuple, Optional, Set, Tuple, Union, cast, TYPE_CHECKING
)

import numpy as np
from sympy.printing.lambdarepr import LambdaPrinter
from sympy.utilities.iterables import iterable, flatten

# Handle type checking to avoid circular import issues
if TYPE_CHECKING:
    from kingdon.multivector import MultiVector

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
    General wrapper to generate code for an operation on multivectors.
    
    This function handles the various types of outputs that codegen functions
    might return and converts them to a consistent CodegenOutput.
    
    Args:
        codegen: The codegen function to use.
        *mvs: Multivectors to operate on.
    
    Returns:
        CodegenOutput with keys and a callable function.
    
    Example:
        >>> gp_result = do_codegen(codegen_gp, vector1, vector2)
        >>> # Use the result to compute the geometric product
        >>> result = gp_result.func(vector1.values(), vector2.values())
    """
    algebra = mvs[0].algebra
    res = codegen(*mvs)
    if isinstance(res, CodegenOutput):
        return res
    if isinstance(res, LambdifyInput):
        funcname = res.funcname
        args = res.args
        dependencies = res.dependencies
        res = res.expr_dict
    else:
        funcname = f'{codegen.__name__}_' + '_x_'.join(f"{mv.type_number}" for mv in mvs)
        args = {arg_name: arg.values() for arg_name, arg in zip(string.ascii_uppercase, mvs)}
        dependencies = None

    # Sort keys in canonical order.
    res = {bin: res[bin] if isinstance(res, dict) else getattr(res, bin)
           for bin in res.keys() if bin in algebra.canon2bin.values()}

    if not algebra.cse and any(isinstance(v, str) for v in res.values()):
        return func_builder(res, *mvs, funcname=funcname)

    keys, exprs = tuple(res.keys()), list(res.values())
    func = lambdify(args, exprs, funcname=funcname, cse=algebra.cse, dependencies=dependencies)
    return CodegenOutput(keys, func)


def do_compile(codegen: Callable, *tapes: Any) -> CodegenOutput:
    """
    Compile a codegen function with already-compiled input "tapes".
    
    This is used for chaining operations at a lower level than do_codegen.
    
    Args:
        codegen: The codegen function to use.
        *tapes: Pre-compiled inputs.
    
    Returns:
        CodegenOutput with keys and a callable function.
    """
    algebra = tapes[0].algebra
    namespace = algebra.numspace
    res = codegen(*tapes)
    funcname = f'{codegen.__name__}_' + '_x_'.join(f"{tape.type_number}" for tape in tapes)
    funcstr = f"def {funcname}({', '.join(t.expr for t in tapes)}):"
    if not isinstance(res, str):
        funcstr += f"    return {res.expr}"
    else:
        funcstr += f"    return ({res},)"
    funclocals = {}
    filename = f'<{funcname}>'
    c = compile(funcstr, filename, 'exec')
    exec(c, namespace, funclocals)
    linecache.cache[filename] = (len(funcstr), None, funcstr.splitlines(True), filename)  # type: ignore
    func = funclocals[funcname]
    return CodegenOutput(
        res.keys() if not isinstance(res, str) else (0,),
        func
    )


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

    def _print_unpacking(self, unpackto: Any, arg: str) -> List[str]:
        """
        Generate code for unpacking arguments into variables.
        
        Args:
            unpackto: Structure to unpack into.
            arg: Argument to unpack.
        
        Returns:
            List of code lines for unpacking.
        """
        def unpack_lhs(lvalues: Any) -> str:
            return f'[{", ".join(unpack_lhs(val) if iterable(val) else val for val in lvalues)}]'
        return [f'{unpack_lhs(unpackto)} = {arg}']


def _recursive_to_string(doprint: Callable, arg: Any) -> str:
    """
    Recursively convert expressions to strings.
    
    Args:
        doprint: Function to convert atomic expressions to strings.
        arg: Expression to convert.
    
    Returns:
        String representation of the expression.
    
    Raises:
        NotImplementedError: If the argument type is not supported.
    """
    if isinstance(arg, str):
        return arg
    elif not arg:
        return str(arg)
    elif iterable(arg):
        if isinstance(arg, list):
            left, right = "[", "]"
        elif isinstance(arg, tuple):
            left, right = "(", ",)"
        else:
            raise NotImplementedError(f"Unhandled type: {type(arg)}, {arg}")
        return f'{left}{", ".join(_recursive_to_string(doprint, e) for e in arg)}{right}'
    else:
        return doprint(arg)


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
    tosympy = lambda x: x.tosympy() if hasattr(x, 'tosympy') else x
    args = {name: [tosympy(v) for v in values] for name, values in args.items()}
    exprs = [tosympy(expr) for expr in exprs]
    if dependencies is not None:
        dependencies = [(tosympy(y), tosympy(x)) for y, x in dependencies]
    
    names = tuple(arg if isinstance(arg, str) else arg.name for arg in args.keys())
    iterable_args = tuple(args.values())
    funcprinter = KingdonPrinter(printer, dummify)
    
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
    
    funcstr = funcprinter.doprint(funcname, iterable_args, names, _exprs, cses=cses)
    namespace = {'builtins': builtins, 'range': range}
    funclocals = {}
    filename = f'<{funcname}>'
    c = compile(funcstr, filename, 'exec')
    exec(c, namespace, funclocals)
    linecache.cache[filename] = (len(funcstr), None, funcstr.splitlines(True), filename)  # type: ignore
    func = funclocals[funcname]
    return func