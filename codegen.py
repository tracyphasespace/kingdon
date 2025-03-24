"""
Code Generation Module for Geometric Algebra
============================================

This module provides code generation utilities for Geometric Algebra operations.
It creates optimized functions for various products (geometric, outer, inner, etc.)
and other operations on multivectors.

This is a wrapper module that imports and re-exports functions from:
- codegen_core.py: Core utilities and infrastructure
- codegen_operations.py: Operation-specific implementations

The module allows for efficient computation of Geometric Algebra operations through
code generation, which dynamically creates optimized functions for specific operations.

Example:
    >>> try:
    try:
    try:
    try:
    try:
    from kingdon.algebra import Algebra
except ImportError:
    # During initial import
    Algebra = object
except ImportError:
    # During initial import
    Algebra = object
except ImportError:
    # During initial import
    Algebra = object
except ImportError:
    # During initial import
    Algebra = object
except ImportError:
    # During initial import
    Algebra = object
    >>> from kingdon.codegen import do_codegen, codegen_gp
    >>> alg = Algebra(3, 0, 1)  # 3D projective geometric algebra
    >>> x = alg.multivector(name="x", grades=(1,))
    >>> y = alg.multivector(name="y", grades=(1,))
    >>> gp_func = do_codegen(codegen_gp, x, y)
    >>> # gp_func can now be used to compute the geometric product of x and y
"""

# Re-export core utilities
from kingdon.codegen_core import (
    mathstr,
    AdditionChains,
    CodegenOutput,
    LambdifyInput,
    Fraction,
    power_supply,
    codegen_product,
    do_codegen,
    do_compile,
    func_builder,
    lambdify,
    KingdonPrinter
)

# Re-export operation implementations
from kingdon.codegen_operations import (
    # Product operations
    codegen_gp,
    codegen_sw,
    codegen_cp,
    codegen_ip,
    codegen_op,
    codegen_div,
    codegen_rp,
    codegen_acp,
    codegen_proj,
    codegen_sp,
    codegen_lc,
    codegen_rc,
    
    # Inverse operations
    codegen_inv,
    codegen_hitzer_inv,
    codegen_shirokov_inv,
    
    # Transformations
    codegen_normsq,
    codegen_reverse,
    codegen_involute,
    codegen_conjugate,
    codegen_involutions,
    codegen_sqrt,
    
    # Duality operations
    codegen_polarity,
    codegen_unpolarity,
    codegen_hodge,
    codegen_unhodge,
    
    # Exponential functions
    codegen_outerexp,
    codegen_outersin,
    codegen_outercos,
    codegen_outertan,
    
    # Basic arithmetic
    codegen_add,
    codegen_sub,
    codegen_neg
)

# Make all imported items available when importing the module
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