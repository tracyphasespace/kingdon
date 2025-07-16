from __future__ import annotations
print("--- DEBUG: Executing with the NEW, CORRECTED codegen.py ---") # <-- CORRECT LOCATION


# Standard Library Imports
import string
import builtins
import keyword
import inspect
import linecache
import operator
import warnings
import math
import sympy
import logging
from collections import Counter, defaultdict, namedtuple
from dataclasses import dataclass
from functools import reduce, cached_property, partial
from itertools import product, combinations, groupby

from sympy.printing.numpy import NumPyPrinter
from sympy import Mul, Add # Make sure Mul and Add are imported
import ast # For literal_eval safely

# Typing Imports
from typing import (
    Any, Callable, Dict, Generator, Iterable, List, Mapping,
    NamedTuple, Optional, Set, Tuple, Union, cast, TYPE_CHECKING
)

# Third-party Imports
import numpy as np
import sympy
# Import specific sympy elements needed
from sympy import Expr, Symbol, sympify, Integer, Add, Mul, Pow, Rational, simplify, expand, factor, collect, Basic # Added Basic for type checks

from sympy.printing.lambdarepr import LambdaPrinter
# Or, if you already have `import sympy.printing`, change it to:
# import sympy.printing.lambdarepr
# and then use sympy.printing.lambdarepr.LambdaPrinter

# Setup logging
log = logging.getLogger(__name__)

# Defer Internal Imports using TYPE_CHECKING
if TYPE_CHECKING:
    from .multivector import MultiVector
    from .algebra import Algebra, AlgebraError # Import AlgebraError if defined there

# Define AlgebraError locally if not imported
# Ensure this matches the definition in operator_dict.py or a central location
try:
    from .operator_dict import AlgebraError
except ImportError:
    class AlgebraError(Exception):
        """Custom exception for errors related to Geometric Algebra operations."""
        pass


# Helper function to safely import MultiVector only when needed at runtime
def _get_mv_class():
    """Dynamically imports and returns the MultiVector class."""
    from .multivector import MultiVector
    return MultiVector

#=============================================================================
# CORE UTILITIES
#=============================================================================

class CodegenOutput(NamedTuple):
    """ Output of a codegen function containing the keys and callable function. """
    # (Docstring unchanged)
    keys_out: Tuple[int, ...]
    func: Callable[..., List[Any]]

def _get_arg_symbols(mv: 'MultiVector', name: str) -> List[Symbol]:
    """ Create SymPy symbols for each component of a multivector. """
    # (Method unchanged)
    algebra = mv.algebra
    return [Symbol(f"{name}_{algebra.bin2canon.get(k, f'k{k}')}") for k in mv.keys()]

def _values_to_symbol_dict(mv: 'MultiVector', symbols: List[Symbol]) -> Dict[int, Symbol]:
    """ Map multivector keys to their corresponding SymPy symbols. """
    # (Method unchanged)
    if len(mv.keys()) != len(symbols):
         base_name = symbols[0].name.split('_')[0] if symbols else 'unknown_mv'
         raise ValueError(f"Mismatch between number of keys ({len(mv.keys())}) and symbols ({len(symbols)}) for MV '{base_name}'.")
    return dict(zip(mv.keys(), symbols))


def codegen_product(x_sym: Dict[int, Expr], y_sym: Dict[int, Expr],
                    algebra: 'Algebra',
                    filter_func: Optional[Callable[[int, int, int], bool]] = None,
                    sign_func: Optional[Callable[[Tuple[int, int]], int]] = None,
                    keyout_func: Callable[[int, int], int] = operator.xor) -> Dict[int, Expr]:
    """ Helper function for symbolic code generation of product-type operations. """
    # (Now accepts Expr for values, not just Symbol)
    sign_func = sign_func or (lambda pair: algebra.signs.get(pair, 0))
    res: Dict[int, Expr] = defaultdict(lambda: Integer(0))
    # Ensure keys are valid before iterating (should be handled by caller)
    x_items = list(x_sym.items())
    y_items = list(y_sym.items())

    for (kx, vx), (ky, vy) in product(x_items, y_items):
        # Ensure vx, vy are sympy expressions
        vx_expr = sympify(vx)
        vy_expr = sympify(vy)

        sign = sign_func((kx, ky))
        if sign != 0:
            key_out = keyout_func(kx, ky)
            if filter_func and not filter_func(kx, ky, key_out):
                continue
            # Create term carefully, avoiding unnecessary evaluation
            term_parts = [Integer(sign)]
            if vx_expr != Integer(1): term_parts.append(vx_expr)
            if vy_expr != Integer(1): term_parts.append(vy_expr)
            term = Mul(*term_parts, evaluate=False)

            # Accumulate terms using sympy Add
            res[key_out] = Add(res[key_out], term, evaluate=False)

    # Simplify final accumulated expressions
    final_res = {}
    for k, v_accumulated in res.items():
         # Use doit() before simp_func for potential cancellations
         v_simplified = algebra.simp_func(v_accumulated.doit())
         # Ensure simplification result is still a SymPy expression
         if not isinstance(v_simplified, Basic): v_simplified = sympify(v_simplified)
         if v_simplified != Integer(0):
              final_res[k] = v_simplified
    return final_res


def do_codegen(codegen_func: Callable[..., Dict[int, Expr]], *mvs: 'MultiVector') -> CodegenOutput:
    """ Generate, compile, and return an optimized function for a GA operation. """
    # (Docstring unchanged)
    if not mvs:
        raise ValueError("At least one multivector must be provided to do_codegen.")

    algebra = mvs[0].algebra
    MultiVector = _get_mv_class()
    op_name = codegen_func.__name__.replace('codegen_', '') # Get operation name
    input_keys_repr = tuple(mv.keys() for mv in mvs) # For error messages

    # --- Generate Symbolic Arguments ---
    arg_names = string.ascii_lowercase[:len(mvs)]
    symbolic_args_list = []
    symbolic_arg_map_list = []
    try:
        for i, mv in enumerate(mvs):
             if not isinstance(mv, MultiVector): # Basic type check
                  raise TypeError(f"Input argument {i} (type: {type(mv).__name__}) is not a MultiVector.")
             if mv.algebra != algebra:
                  raise ValueError(f"Input MVs must belong to the same algebra (expected {algebra}, got {mv.algebra} for arg {i}).")
             arg_symbols = _get_arg_symbols(mv, arg_names[i])
             symbolic_args_list.append(arg_symbols)
             symbolic_arg_map_list.append(_values_to_symbol_dict(mv, arg_symbols))
    except (ValueError, TypeError) as e:
         raise AlgebraError(f"Error preparing symbolic arguments for '{op_name}' with keys {input_keys_repr}: {e}") from e

    # --- Call Symbolic Generator ---
    sig = inspect.signature(codegen_func)
    codegen_kwargs = {}
    if 'algebra' in sig.parameters: codegen_kwargs['algebra'] = algebra
    try:
        symbolic_result_dict = codegen_func(*symbolic_arg_map_list, **codegen_kwargs)
        if symbolic_result_dict is NotImplemented:
             raise NotImplementedError(f"Symbolic codegen for '{op_name}' is not implemented for input keys {input_keys_repr}.")
    except NotImplementedError: raise # Propagate explicitly
    except (ValueError, TypeError, ZeroDivisionError) as e_codegen:
        # Catch specific errors expected from codegen_* funcs
        log.warning(f"Symbolic generation for '{op_name}' failed for keys {input_keys_repr}: {e_codegen}")
        # Re-raise as AlgebraError to indicate failure in this specific operation/input
        raise AlgebraError(f"Symbolic generation for '{op_name}' failed for keys {input_keys_repr}: {e_codegen}") from e_codegen
    except Exception as e_codegen:
        # Catch unexpected errors during symbolic generation
        log.error(f"Unexpected error during symbolic generation in {op_name} for keys {input_keys_repr}: {e_codegen}", exc_info=True)
        raise AlgebraError(f"Unexpected error during symbolic generation for '{op_name}' with keys {input_keys_repr}: {e_codegen}") from e_codegen

    # --- Process Symbolic Result ---
    try:
        if not isinstance(symbolic_result_dict, dict):
             if symbolic_result_dict == 0: symbolic_result_dict = {0: Integer(0)}
             elif isinstance(symbolic_result_dict, (int, float, complex, Expr)):
                  symbolic_result_dict = {0: sympify(symbolic_result_dict)}
             else: raise TypeError(f"{codegen_func.__name__} returned unexpected type {type(symbolic_result_dict).__name__}, expected dict.")

        # Filter zero results after simplification
        final_symbolic_result_dict = {
            k: v_simp for k, v in symbolic_result_dict.items()
            if not (v_simp := algebra.simp_func(v)).is_zero # Use simp_func and check is_zero
        }
    except Exception as e_process:
         raise AlgebraError(f"Error processing symbolic result from '{op_name}' for keys {input_keys_repr}: {e_process}") from e_process

    # --- Determine Output Keys and Expressions ---
    if not final_symbolic_result_dict:
        keys_out = (0,)
        exprs_out = [Integer(0)] # Use SymPy zero
    else:
        # Sorting logic remains the same
        def sort_key(k):
             grade = bin(k).count('1')
             try: grade_indices = algebra.indices_for_grades.get((grade,), []); return (grade, grade_indices.index(k))
             except (ValueError, IndexError, KeyError, AttributeError): return (grade, k)
        keys_out = tuple(sorted(final_symbolic_result_dict.keys(), key=sort_key))
        exprs_out = [final_symbolic_result_dict[k] for k in keys_out]

    # --- Compile function using lambdify ---
    # Create a descriptive function name
    keys_repr = '_'.join(f"{''.join(map(str,k))}" for k in input_keys_repr) # Compact key representation
    funcname = f'{op_name}_{keys_repr}'
    all_input_symbols = [symbol for sublist in symbolic_args_list for symbol in sublist]

    try:
        compiled_func = lambdify(
            all_input_symbols, exprs_out, funcname=funcname, cse=algebra.cse,
        )
        log.debug(f"Successfully compiled {funcname} for input keys: {input_keys_repr}")
    except Exception as e_compile:
        log.error(f"Lambdify compilation error for {funcname} with symbolic expressions {exprs_out}: {e_compile}", exc_info=True)
        # Provide a non-functional placeholder that raises an error
        def error_func(*args, **kwargs):
            raise AlgebraError(f"Codegen compilation failed for {funcname} (operation: {op_name}, input keys: {input_keys_repr}). Lambdify error: {e_compile}")
        compiled_func = error_func
        keys_out = keys_out if keys_out else (0,) # Ensure keys_out is valid

    # In codegen.py, inside the `do_codegen` function

    # ... (code before the wrapper function) ...

    num_expected_args = len(mvs) # Number of original symbolic multivectors

# In codegen.py, inside the `do_codegen` function
# THIS IS THE ORIGINAL, CORRECT WRAPPER

# In codegen.py, inside the `do_codegen` function
# THIS IS THE CORRECT, SIMPLIFIED WRAPPER

    def wrapper_func(*runtime_args):
        """
        Wrapper to execute the compiled function and handle runtime errors.
        It receives a flat tuple of all coefficient values.
        """
        # The compiled function expects a flat sequence of arguments.
        if len(runtime_args) != len(all_input_symbols):
             raise TypeError(f"{funcname} expected {len(all_input_symbols)} arguments, got {len(runtime_args)}")
        
        try:
            # Execute the JIT-compiled function
            result = compiled_func(*runtime_args)
            # Ensure the output is always a list, as expected by the caller
            return result if isinstance(result, list) else [result]
        except Exception as e_runtime:
             # Add context to any error that occurs inside the compiled code
             raise AlgebraError(f"Runtime error in compiled function '{funcname}': {e_runtime}") from e_runtime

    wrapper_func.__name__ = f"{funcname}_runtime_wrapper"
    return CodegenOutput(keys_out, wrapper_func)






def lambdify(args: List[Symbol], exprs: List[Expr], funcname: str,
             modules: Optional[List[Union[str, Dict]]] = None,
             printer: Any = None,
             dummify: bool = False, cse: bool = False) -> Callable:
    """ Turn symbolic expressions into a Python function using SymPy's lambdify. """

    if printer is None:
         printer = NumPyPrinter()

    # THE DEFINITIVE FIX IS HERE:
    # We are explicitly mapping the symbolic SymPy operations `Mul` and `Add`
    # to the desired NumPy functions `np.multiply` and `np.add`.
    # This leaves no ambiguity for the `*` and `+` operators.
    if modules is None:
        modules = [
            {'Mul': np.multiply, 'Add': np.add}, # Explicit mapping
            'numpy'                             # General numpy namespace
        ]

    try:
         # The rest of this function remains the same as before.
         valid_exprs = []
         for i, expr in enumerate(exprs):
              if isinstance(expr, Basic): valid_exprs.append(expr)
              else:
                   try: valid_exprs.append(sympify(expr))
                   except (TypeError, ValueError) as e_sympify:
                        raise TypeError(f"Expression at index {i} ('{expr}') cannot be converted to SymPy expression for lambdify in '{funcname}'. Error: {e_sympify}") from e_sympify

         compiled_func = sympy.lambdify(args, valid_exprs, modules=modules, printer=printer, dummify=dummify, cse=cse)

         # Add source code to linecache for better tracebacks
         filename = f'<lambdify-{funcname}>'
         compiled_func.__code__ = compiled_func.__code__.replace(co_filename=filename)
         func_sig = f"def {funcname}({', '.join(map(str, args))}):"
         try:
              if isinstance(valid_exprs, (list, tuple)):
                  if len(valid_exprs) == 1: return_stmt = f"    return {printer.doprint(valid_exprs[0])}"
                  else: return_stmt = f"    return [{', '.join(printer.doprint(e) for e in valid_exprs)}]"
              else: return_stmt = f"    return {printer.doprint(valid_exprs)}"
              fake_source = f"{func_sig}\n{return_stmt}"
              linecache.cache[filename] = (len(fake_source), None, fake_source.splitlines(True), filename)
         except Exception as e_print:
              log.warning(f"Could not generate source representation for {funcname}: {e_print}")
              linecache.cache[filename] = (1, None, [f"# Error generating source for {funcname}"], filename)

         return compiled_func
    except Exception as e:
         log.error(f"SymPy lambdify failed for '{funcname}' with args {args} and expressions {exprs}: {e}", exc_info=True)
         raise AlgebraError(f"Failed to compile function '{funcname}' using lambdify. Error: {e}") from e




#=============================================================================
# OPERATION-SPECIFIC SYMBOLIC GENERATORS
#=============================================================================
# (Focus on raising specific ValueErrors for logical failures,
#  and propagating NotImplementedError for missing symbolic logic)

# --- Helper for Symbolic Dictionary Multiplication ---
def _symbolic_dict_mul(res_dict: Dict[int, Expr], factor: Expr, algebra: 'Algebra') -> Dict[int, Expr]:
    # (Method unchanged)
    out_dict = {}; factor_simp = algebra.simp_func(factor)
    if factor_simp == Integer(0): return {0: Integer(0)}
    for k, v in res_dict.items():
        term = Mul(factor_simp, v, evaluate=False); term_simplified = algebra.simp_func(term.doit())
        if not isinstance(term_simplified, Basic): term_simplified = sympify(term_simplified) # Ensure sympy type
        if term_simplified != Integer(0): out_dict[k] = term_simplified
    return out_dict if out_dict else {0: Integer(0)}

# --- Product Type Operations ---
def codegen_gp(x_sym: Dict[int, Expr], y_sym: Dict[int, Expr], algebra: 'Algebra') -> Dict[int, Expr]:
    """ Generate symbolic dictionary for the geometric product (x * y). """
    return codegen_product(x_sym, y_sym, algebra=algebra)

def _outer_product_filter(kx: int, ky: int, k_out: int) -> bool:
    # (Function unchanged)
    grade_kx = bin(kx).count('1'); grade_ky = bin(ky).count('1'); grade_k_out = bin(k_out).count('1')
    return grade_k_out == grade_kx + grade_ky # Outer product grade condition (Removed kx&ky==0 which is implicit)

def codegen_op(x_sym: Dict[int, Expr], y_sym: Dict[int, Expr], algebra: 'Algebra') -> Dict[int, Expr]:
    """ Generate symbolic dictionary for the outer product (x ^ y). """
    return codegen_product(x_sym, y_sym, algebra=algebra, filter_func=_outer_product_filter)

def codegen_ip(x_sym: Dict[int, Expr], y_sym: Dict[int, Expr], algebra: 'Algebra', diff_func: Callable[[int], int] = abs) -> Dict[int, Expr]:
    """ Generate symbolic dictionary for an inner product (x | y). """
    filter_func = lambda kx, ky, k_out: (bin(k_out).count('1') == diff_func(bin(kx).count('1') - bin(ky).count('1')))
    return codegen_product(x_sym, y_sym, algebra=algebra, filter_func=filter_func)

def codegen_lc(x_sym: Dict[int, Expr], y_sym: Dict[int, Expr], algebra: 'Algebra') -> Dict[int, Expr]:
    """ Generate symbolic dictionary for the left contraction (x << y). """
    return codegen_ip(x_sym, y_sym, algebra=algebra, diff_func=lambda k: -k)

def codegen_rc(x_sym: Dict[int, Expr], y_sym: Dict[int, Expr], algebra: 'Algebra') -> Dict[int, Expr]:
    """ Generate symbolic dictionary for the right contraction (x >> y). """
    return codegen_ip(x_sym, y_sym, algebra=algebra, diff_func=lambda k: k)

def codegen_sp(x_sym: Dict[int, Expr], y_sym: Dict[int, Expr], algebra: 'Algebra') -> Dict[int, Expr]:
    """ Generate symbolic dictionary for the scalar product (x . y). """
    return codegen_ip(x_sym, y_sym, algebra=algebra, diff_func=lambda k: 0)

# --- Commutators ---
def codegen_cp(x_sym: Dict[int, Expr], y_sym: Dict[int, Expr], algebra: 'Algebra') -> Dict[int, Expr]:
    """ Generate symbolic dictionary for the commutator product 0.5 * (x*y - y*x). """
    try:
        xy = codegen_gp(x_sym, y_sym, algebra); yx = codegen_gp(y_sym, x_sym, algebra)
    except NotImplementedError as nie: raise NotImplementedError("Symbolic CP requires symbolic GP.") from nie
    result = defaultdict(lambda: Integer(0)); all_keys = set(xy.keys()) | set(yx.keys()); half = Rational(1, 2)
    for k in all_keys:
        term = half * (xy.get(k, Integer(0)) - yx.get(k, Integer(0)))
        term_simplified = algebra.simp_func(term)
        if not isinstance(term_simplified, Basic): term_simplified = sympify(term_simplified)
        if term_simplified != Integer(0): result[k] = term_simplified
    return dict(result)

def codegen_acp(x_sym: Dict[int, Expr], y_sym: Dict[int, Expr], algebra: 'Algebra') -> Dict[int, Expr]:
    """ Generate symbolic dictionary for the anti-commutator product 0.5 * (x*y + y*x). """
    try:
        xy = codegen_gp(x_sym, y_sym, algebra); yx = codegen_gp(y_sym, x_sym, algebra)
    except NotImplementedError as nie: raise NotImplementedError("Symbolic ACP requires symbolic GP.") from nie
    result = defaultdict(lambda: Integer(0)); all_keys = set(xy.keys()) | set(yx.keys()); half = Rational(1, 2)
    for k in all_keys:
        term = half * (xy.get(k, Integer(0)) + yx.get(k, Integer(0)))
        term_simplified = algebra.simp_func(term)
        if not isinstance(term_simplified, Basic): term_simplified = sympify(term_simplified)
        if term_simplified != Integer(0): result[k] = term_simplified
    return dict(result)

# --- Other Products ---
def codegen_rp(x_sym: Dict[int, Expr], y_sym: Dict[int, Expr], algebra: 'Algebra') -> Dict[int, Expr]:
    """ Generate symbolic dictionary for the regressive product (x & y) using Hodge dual. """
    log.debug("Generating symbolic regressive product using Hodge dual.")
    try:
        dual_x = codegen_hodge(x_sym, algebra); dual_y = codegen_hodge(y_sym, algebra)
        op_duals = codegen_op(dual_x, dual_y, algebra)
        result = codegen_unhodge(op_duals, algebra)
        return result
    except NotImplementedError as nie: raise NotImplementedError(f"Symbolic RP requires working dependencies (Hodge, OP, Unhodge): {nie}") from nie
    except (ValueError, TypeError) as e: raise AlgebraError(f"Error during symbolic RP generation: {e}") from e

def codegen_sw(x_sym: Dict[int, Expr], y_sym: Dict[int, Expr], algebra: 'Algebra') -> Dict[int, Expr]:
    r""" Generate symbolic dictionary for the sandwich product (x * y * ~x). """
    log.warning("Symbolic codegen_sw requires robust symbolic dependencies (GP, REV).")
    try:
         xy = codegen_gp(x_sym, y_sym, algebra)
         x_rev = codegen_reverse(x_sym, algebra)
         result = codegen_gp(xy, x_rev, algebra)
         return result
    except NotImplementedError as nie:
         raise NotImplementedError(f"Symbolic SW requires symbolic GP and REV: {nie}") from nie
    except (ValueError, TypeError) as e:
         raise AlgebraError(f"Error during symbolic SW generation: {e}") from e


def codegen_proj(x_sym: Dict[int, Expr], y_sym: Dict[int, Expr], algebra: 'Algebra') -> Dict[int, Expr]:
    r""" Generate symbolic dictionary for projection of x onto y: (x | y) * y.inv(). """
    log.debug("Generating symbolic projection using IP and Inverse.")
    try:
        ip_xy = codegen_ip(x_sym, y_sym, algebra=algebra)
        # If IP is zero, projection is zero
        if not ip_xy or all(v.is_zero for v in ip_xy.values()):
            return {0: Integer(0)}
        # Try to compute inverse of y
        inv_y = codegen_inv(y_sym, algebra=algebra) # Can raise ValueError/NotImplementedError

        result = codegen_gp(ip_xy, inv_y, algebra=algebra) # GP now handles Expr inputs
        return result
    except NotImplementedError as nie: raise NotImplementedError(f"Symbolic projection requires symbolic IP, INV, GP: {nie}") from nie
    except ValueError as ve: raise AlgebraError(f"Cannot compute symbolic projection: failed to invert target y. {ve}") from ve
    except TypeError as te: raise AlgebraError(f"Type error during symbolic projection: {te}") from te


# --- Basic Arithmetic ---
def codegen_add(x_sym: Dict[int, Expr], y_sym: Dict[int, Expr], algebra: 'Algebra') -> Dict[int, Expr]:
    """ Generate symbolic dictionary for addition (x + y). """
    res = defaultdict(lambda: Integer(0)); all_keys = set(x_sym.keys()) | set(y_sym.keys())
    for k in all_keys:
        term = Add(x_sym.get(k, Integer(0)), y_sym.get(k, Integer(0)), evaluate=False)
        res[k] = term # Accumulate symbolically
    final_res = {}
    for k, v_accumulated in res.items():
         v_simplified = algebra.simp_func(v_accumulated.doit())
         if not isinstance(v_simplified, Basic): v_simplified = sympify(v_simplified)
         if v_simplified != Integer(0): final_res[k] = v_simplified
    return final_res

def codegen_sub(x_sym: Dict[int, Expr], y_sym: Dict[int, Expr], algebra: 'Algebra') -> Dict[int, Expr]:
    """ Generate symbolic dictionary for subtraction (x - y). """
    res = defaultdict(lambda: Integer(0)); all_keys = set(x_sym.keys()) | set(y_sym.keys())
    for k in all_keys:
        term = Add(x_sym.get(k, Integer(0)), -sympify(y_sym.get(k, Integer(0))), evaluate=False)
        res[k] = term
    final_res = {}
    for k, v_accumulated in res.items():
         v_simplified = algebra.simp_func(v_accumulated.doit())
         if not isinstance(v_simplified, Basic): v_simplified = sympify(v_simplified)
         if v_simplified != Integer(0): final_res[k] = v_simplified
    return final_res

def codegen_neg(x_sym: Dict[int, Expr], algebra: 'Algebra') -> Dict[int, Expr]:
    """ Generate symbolic dictionary for negation (-x). """
    res = {};
    for k, v in x_sym.items():
        term = -sympify(v); term_simplified = algebra.simp_func(term)
        if not isinstance(term_simplified, Basic): term_simplified = sympify(term_simplified)
        if term_simplified != Integer(0): res[k] = term_simplified
    return res

# --- Unary Transformations ---
def codegen_reverse(x_sym: Dict[int, Expr], algebra: 'Algebra') -> Dict[int, Expr]:
    """ Generate symbolic dictionary for reversion (~x). """
    res = {};
    for k, v in x_sym.items():
        grade = bin(k).count('1'); sign = -1 if grade % 4 in (2, 3) else 1
        term = Mul(sign, sympify(v), evaluate=False); term_simplified = algebra.simp_func(term.doit())
        if not isinstance(term_simplified, Basic): term_simplified = sympify(term_simplified)
        if term_simplified != Integer(0): res[k] = term_simplified
    return res

def codegen_involute(x_sym: Dict[int, Expr], algebra: 'Algebra') -> Dict[int, Expr]:
    """ Generate symbolic dictionary for grade involution (x*). """
    res = {};
    for k, v in x_sym.items():
        grade = bin(k).count('1'); sign = -1 if grade % 2 != 0 else 1
        term = Mul(sign, sympify(v), evaluate=False); term_simplified = algebra.simp_func(term.doit())
        if not isinstance(term_simplified, Basic): term_simplified = sympify(term_simplified)
        if term_simplified != Integer(0): res[k] = term_simplified
    return res

def codegen_conjugate(x_sym: Dict[int, Expr], algebra: 'Algebra') -> Dict[int, Expr]:
    """ Generate symbolic dictionary for Clifford conjugation (x†). """
    res = {};
    for k, v in x_sym.items():
        grade = bin(k).count('1'); sign = -1 if grade % 4 in (1, 2) else 1
        term = Mul(sign, sympify(v), evaluate=False); term_simplified = algebra.simp_func(term.doit())
        if not isinstance(term_simplified, Basic): term_simplified = sympify(term_simplified)
        if term_simplified != Integer(0): res[k] = term_simplified
    return res

# --- Inverse, Division, Norm ---
def codegen_inv(y_sym: Dict[int, Expr], algebra: 'Algebra') -> Dict[int, Expr]:
     """ Generate symbolic dictionary for the inverse y.inv() = ~y / (y * ~y). """
     log.debug("Generating symbolic inverse using ~y / (y*~y).")
     try:
         y_rev_sym = codegen_reverse(y_sym, algebra)
         y_ytilde_sym = codegen_gp(y_sym, y_rev_sym, algebra)

         if not y_ytilde_sym: # Check if product is identically zero
             raise ValueError("Cannot invert multivector: y*~y is symbolically zero.")

         # Extract and simplify scalar part
         norm_sq_expr = algebra.simp_func(y_ytilde_sym.get(0, Integer(0)).doit())
         if not isinstance(norm_sq_expr, Basic): norm_sq_expr = sympify(norm_sq_expr)

         # Check if simplified scalar part is zero
         if norm_sq_expr.is_zero:
             raise ValueError("Cannot invert multivector: scalar part of y*~y simplifies to zero.")

         # Check for non-scalar terms after simplification
         has_non_scalar_terms = False
         for k, v in y_ytilde_sym.items():
              if k != 0:
                  v_simplified = algebra.simp_func(v.doit())
                  if not isinstance(v_simplified, Basic): v_simplified = sympify(v_simplified)
                  if not v_simplified.is_zero:
                       has_non_scalar_terms = True
                       log.warning(f"Symbolic inverse: y*~y has non-scalar component for key {k}: {v_simplified}")
                       break # Found one non-scalar term

         if has_non_scalar_terms:
              raise ValueError("Symbolic inverse failed: y*~y has non-scalar components.")

         # Check if the scalar part is truly scalar (independent of basis elements)
         # Note: This check might be complex for general expressions. is_constant() is often sufficient.
         if not norm_sq_expr.is_constant():
              warnings.warn(f"Symbolic inverse: y*~y simplified to non-constant expression '{norm_sq_expr}'. Inversion might be complex.", stacklevel=3)
              # Proceed, but result might be complicated

         # Perform division: ~y * (1 / norm_sq_expr)
         try: inv_norm_sq = Pow(norm_sq_expr, -1, evaluate=False) # Use Pow for symbolic inverse
         except Exception as e_pow: raise AlgebraError(f"Error creating symbolic power for inverse norm^2: {e_pow}") from e_pow

         result_dict = _symbolic_dict_mul(y_rev_sym, inv_norm_sq, algebra)
         return result_dict

     except NotImplementedError as nie: raise NotImplementedError(f"Symbolic INV requires symbolic REV and GP: {nie}") from nie
     except ValueError as ve: raise # Re-raise specific value errors
     except Exception as e: raise AlgebraError(f"Unexpected error during symbolic inverse: {e}") from e


def codegen_div(x_sym: Dict[int, Expr], y_sym: Dict[int, Expr], algebra: 'Algebra') -> Dict[int, Expr]:
    """ Generates symbolic division x / y = x * y.inv(). """
    log.debug("Generating symbolic division using x * y.inv().")
    try:
        inv_y = codegen_inv(y_sym, algebra=algebra) # Can raise ValueError/NotImplementedError/AlgebraError
        result = codegen_gp(x_sym, inv_y, algebra=algebra)
        return result
    except NotImplementedError as nie: raise NotImplementedError(f"Symbolic DIV requires symbolic GP and INV: {nie}") from nie
    except (ValueError, AlgebraError) as ve: raise AlgebraError(f"Cannot compute symbolic division: Denominator 'y' failed inversion. {ve}") from ve
    except Exception as e: raise AlgebraError(f"Unexpected error during symbolic division: {e}") from e


# Inside kingdon/codegen.py

# ... (other codegen functions) ...

def codegen_normsq(x_sym: Dict[int, Expr], algebra: 'Algebra') -> Dict[int, Expr]:
    """
    Generate symbolic dictionary for squared norm |x|^2 = <x * ~x>_0.
    Corrected to only return the scalar part.
    """
    log.debug("Generating symbolic norm squared using <x * ~x>_0.")
    try:
        # Calculate x * ~x symbolically
        x_rev_sym = codegen_reverse(x_sym, algebra)
        result_dict_full = codegen_gp(x_sym, x_rev_sym, algebra)

        # Extract and simplify only the scalar (grade 0) part
        scalar_part_expr = result_dict_full.get(0, Integer(0)) # Get key 0, default to SymPy 0
        scalar_part_simplified = algebra.simp_func(scalar_part_expr.doit()) # Simplify

        # Ensure the result is a SymPy Basic type
        if not isinstance(scalar_part_simplified, Basic):
            scalar_part_simplified = sympify(scalar_part_simplified)

        # Return a dictionary containing only the scalar part (key 0)
        # Return {0: 0} if the scalar part simplified to zero
        if scalar_part_simplified.is_zero:
            return {0: Integer(0)}
        else:
            return {0: scalar_part_simplified}

    except NotImplementedError as nie:
        # Propagate if dependencies (reverse, gp) are not implemented
        raise NotImplementedError(f"Symbolic normsq requires symbolic REV and GP: {nie}") from nie
    except Exception as e:
        # Catch other unexpected errors
        raise AlgebraError(f"Unexpected error during symbolic norm squared calculation: {e}") from e

# ... (rest of codegen.py) ...


# --- Duals ---
def _numeric_mv_to_symbolic_dict(mv: 'MultiVector', algebra: 'Algebra') -> Dict[int, Expr]:
    """ Converts a numeric MultiVector to a symbolic dict with constant values. """
    # (Ensure values are sympified)
    return {k: sympify(v) for k, v in mv.items()}

def codegen_polarity(x_sym: Dict[int, Expr], algebra: 'Algebra', undual: bool = False) -> Dict[int, Expr]:
    """ Generate symbolic dictionary for polarity (un)dual: x * PSS.inv() or x * PSS. """
    log.debug(f"Generating symbolic polarity dual (undual={undual}).")
    try:
        if undual:
            pss_mv = algebra.pss; pss_target = pss_mv
            if pss_mv is None: raise ValueError("Algebra lacks pseudoscalar (pss) for polarity undual.")
        else:
            try: pss_inv_mv = algebra.pss_inv; pss_target = pss_inv_mv
            except (AttributeError, ZeroDivisionError, ValueError) as e_pss_inv:
                raise ValueError(f"Cannot compute polarity dual: Failed to get/compute PSS inverse. {e_pss_inv}") from e_pss_inv
            if pss_inv_mv is None: raise ValueError("Could not obtain PSS inverse for polarity dual.")

        pss_sym = _numeric_mv_to_symbolic_dict(pss_target, algebra)
        result = codegen_gp(x_sym, pss_sym, algebra)
        return result

    except NotImplementedError as nie: raise NotImplementedError(f"Symbolic polarity requires symbolic GP: {nie}") from nie
    except ValueError as ve: raise # Re-raise value errors (null pss etc)
    except Exception as e: raise AlgebraError(f"Unexpected error during symbolic polarity: {e}") from e


def codegen_unpolarity(x_sym: Dict[int, Expr], algebra: 'Algebra') -> Dict[int, Expr]:
    """ Generate symbolic dictionary for unpolarity (undual using PSS). """
    return codegen_polarity(x_sym, algebra, undual=True)

def codegen_hodge(x_sym: Dict[int, Expr], algebra: 'Algebra', undual: bool = False) -> Dict[int, Expr]:
    """ Generate symbolic dictionary for Hodge (un)dual. """
    # (Method largely unchanged, ensure result values are simplified)
    key_pss = (1 << algebra.d) - 1; res = {}
    sign_func = lambda pair: algebra.signs.get(pair, 0)
    for k_in, v_sym in x_sym.items():
        k_out = key_pss ^ k_in
        sign_pair = (k_in, k_out) if not undual else (k_out, k_in)
        sign = sign_func(sign_pair)
        # Ensure v_sym is sympy type before multiplying
        term = Mul(sign, sympify(v_sym), evaluate=False)
        term_simplified = algebra.simp_func(term.doit())
        if not isinstance(term_simplified, Basic): term_simplified = sympify(term_simplified)
        if term_simplified != Integer(0): res[k_out] = term_simplified
    return res if res else {0: Integer(0)} # Return zero scalar if result is empty

def codegen_unhodge(x_sym: Dict[int, Expr], algebra: 'Algebra') -> Dict[int, Expr]:
    """ Generate symbolic dictionary for Hodge undual. """
    return codegen_hodge(x_sym, algebra, undual=True)

# --- Other Unary (Sqrt, Outer Exp/Sin/Cos/Tan) ---
# (Keep NotImplementedError with clear messages for these)
def codegen_sqrt(x_sym: Dict[int, Expr], algebra: 'Algebra') -> Dict[int, Expr]:
     log.warning("Symbolic codegen_sqrt is complex and generally not supported.")
     raise NotImplementedError("Symbolic square root (codegen_sqrt) is not implemented.")
def codegen_outerexp(x_sym: Dict[int, Expr], algebra: 'Algebra') -> Dict[int, Expr]:
     log.warning("Symbolic outer exponential is complex and not implemented.")
     raise NotImplementedError("Symbolic outer exponential (codegen_outerexp) is not implemented.")
def codegen_outersin(x_sym: Dict[int, Expr], algebra: 'Algebra') -> Dict[int, Expr]:
     log.warning("Symbolic outer sine is complex and not implemented.")
     raise NotImplementedError("Symbolic outer sine (codegen_outersin) is not implemented.")
def codegen_outercos(x_sym: Dict[int, Expr], algebra: 'Algebra') -> Dict[int, Expr]:
     log.warning("Symbolic outer cosine is complex and not implemented.")
     raise NotImplementedError("Symbolic outer cosine (codegen_outercos) is not implemented.")
def codegen_outertan(x_sym: Dict[int, Expr], algebra: 'Algebra') -> Dict[int, Expr]:
     log.warning("Symbolic outer tangent is complex and not implemented.")
     raise NotImplementedError("Symbolic outer tangent (codegen_outertan) is not implemented.")

# Keep __all__ updated
__all__ = [
    'CodegenOutput', 'do_codegen', 'lambdify', 'AlgebraError', # Added AlgebraError
    'codegen_gp', 'codegen_sw', 'codegen_cp', 'codegen_ip', 'codegen_op',
    'codegen_div', 'codegen_rp', 'codegen_acp', 'codegen_proj', 'codegen_sp',
    'codegen_lc', 'codegen_rc',
    'codegen_inv',
    'codegen_normsq',
    'codegen_reverse', 'codegen_involute', 'codegen_conjugate',
    'codegen_sqrt',
    'codegen_polarity', 'codegen_unpolarity', 'codegen_hodge', 'codegen_unhodge',
    'codegen_outerexp', 'codegen_outersin', 'codegen_outercos', 'codegen_outertan',
    'codegen_add', 'codegen_sub', 'codegen_neg'
]