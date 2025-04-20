# -*- coding: utf-8 -*-\
"""
Unified Algebra Module for Geometric Algebra
<<<<<<< HEAD
(Docstring unchanged)
=======
============================================

This module provides a unified Algebra class with core functionality
and factory methods for Geometric Algebra. It consolidates the previous
algebra_base.py and algebra_factory.py into a single module.

The Algebra class represents a Geometric (Clifford) Algebra with signature
(p,q,r) where:
- p: Number of dimensions that square to +1
- q: Number of dimensions that square to -1
- r: Number of dimensions that square to 0 (null/degenerate dimensions)

Key functionality:
- Creation of common algebras (PGA, CGA, VGA, etc.)
- Creation of multivectors (scalars, vectors, bivectors, etc.)
- Geometric operations (products, inversions, duals, etc.)
- Custom operator registry
- Blade dictionary access

Example:
    >>> from kingdon.algebra import Algebra
    >>> # Create a 3D Euclidean geometric algebra
    >>> alg = Algebra(p=3, q=0, r=0)
    >>> # Create a vector
    >>> v = alg.vector([1, 2, 3])
    >>> # Create a bivector
    >>> B = alg.bivector([1, 0, 3])
    >>> # Compute the geometric product
    >>> result = v * B
>>>>>>> 5d1776bf2c10cd778122088d9664260375c22a81
"""

from __future__ import annotations

<<<<<<< HEAD
# Standard Library Imports
import operator
import re
import warnings
import logging
import sympy
import math # Added for factory methods
from itertools import product, combinations, groupby
from functools import partial, reduce, lru_cache, wraps # wraps needed for register
from collections import Counter, defaultdict
from dataclasses import dataclass, field
import inspect # Added for register placeholder

from sympy import sympify
# Or import the whole module:
# import sympy
# and then use sympy.sympify(...)

# Typing Imports
from typing import List, Dict, Optional, Any, Tuple, Callable, Union, Set, Type, TypeVar, cast, TYPE_CHECKING, ClassVar, Sequence # Added Sequence

# Third-party Imports
import numpy as np
import sympy

# Setup logging
log = logging.getLogger(__name__)

# Use TYPE_CHECKING for imports only needed for type hints to avoid cycles
if TYPE_CHECKING:
    from kingdon.multivector import MultiVector
    from kingdon.codegen import CodegenOutput
    from kingdon.operator_dict import OperatorDict, UnaryOperatorDict, BladeDict, Registry, AlgebraError
    # Conditional import for GraphWidget if needed for type hints
    # from kingdon.graph import GraphWidget

# Type variable for operator dictionaries
T = TypeVar('T', bound='OperatorDict')
=======
import operator
import re
import warnings
from itertools import product
from functools import partial, reduce, lru_cache
from collections import Counter, defaultdict
from dataclasses import dataclass, field

from typing import List, Dict, Optional, Any, Tuple, Callable, Union, Set, Type, TypeVar, cast

import numpy as np
import sympy

# Import MultiVector for type annotations
# Use TYPE_CHECKING to avoid circular import issues during runtime
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from kingdon.multivector import MultiVector

# Import code generation functions and supporting classes
from kingdon.codegen import (
    codegen_gp, codegen_sw, codegen_cp, codegen_ip, codegen_op, codegen_div,
    codegen_rp, codegen_acp, codegen_proj, codegen_sp, codegen_lc, codegen_inv,
    codegen_rc, codegen_normsq, codegen_add, codegen_sub, codegen_neg, codegen_reverse,
    codegen_involute, codegen_conjugate, codegen_sqrt, codegen_outerexp, codegen_outersin,
    codegen_outercos, codegen_outertan, codegen_polarity, codegen_unpolarity,
    codegen_hodge, codegen_unhodge
)
from kingdon.operator_dict import OperatorDict, UnaryOperatorDict, BladeDict, Registry, AlgebraError

# Type variable for operator dictionaries
T = TypeVar('T', bound=OperatorDict)
>>>>>>> 5d1776bf2c10cd778122088d9664260375c22a81

# Field initializer for operator dictionaries
operation_field = partial(field, default_factory=dict, init=False, repr=False, compare=False)


@dataclass
class Algebra:
<<<<<<< HEAD
    """ Comprehensive Algebra class for Geometric (Clifford) Algebra. """
=======
    """
    Comprehensive Algebra class for Geometric (Clifford) Algebra.

    This class represents a Geometric Algebra with signature (p,q,r) and provides
    methods for creating multivectors and performing geometric operations.
    
    Attributes:
        p: Number of dimensions that square to +1
        q: Number of dimensions that square to -1
        r: Number of dimensions that square to 0 (null/degenerate dimensions)
        signature: Optional explicit signature as list of 1, -1, and 0 values
        start_index: Starting index for basis vectors (default is 1, or 0 for PGA)
        basis: Optional list of custom basis names
        cse: Whether to use common subexpression elimination
        graded: Whether to use graded algebra representation
        wrapper: Optional function wrapper for code generation
        simp_func: Function for symbolic simplification
        codegen_symbolcls: Symbol class for code generation
        pretty_blade: Pretty string for basis blades
    """
>>>>>>> 5d1776bf2c10cd778122088d9664260375c22a81
    p: int = 0
    q: int = 0
    r: int = 0
    signature: Optional[List[int]] = None
    start_index: Optional[int] = None
<<<<<<< HEAD
    basis_names: Optional[List[str]] = None

    # Internal attributes
    d: int = field(init=False)
    cse: bool = field(default=True, repr=False)
    graded: bool = field(default=False, repr=False)
    simp_func: Callable[[Any], Any] = field(default=lambda x: sympy.simplify(sympy.expand(x)) if isinstance(x, sympy.Expr) else x, repr=False)
    wrapper: Optional[Callable] = field(default=None, repr=False)
    codegen_symbolcls: Any = field(default=sympy.Symbol, repr=False)
    pretty_blade: str = field(default='𝐞', repr=False)

    # Computed mappings and properties - Initialize as empty/defaults
    canon2bin: Dict[str, int] = field(init=False, repr=False, compare=False, default_factory=dict)
    bin2canon: Dict[int, str] = field(init=False, repr=False, compare=False, default_factory=dict)
    _bin2canon_prettystr: Dict[int, str] = field(init=False, repr=False, compare=False, default_factory=dict)
    indices_for_grades: Dict[Tuple[int, ...], List[int]] = field(init=False, repr=False, default_factory=dict)
    signs: Dict[Tuple[int, int], int] = field(init=False, repr=False, compare=False, default_factory=dict)
    blades: 'BladeDict' = field(init=False, repr=False, compare=False)
    pss: 'MultiVector' = field(init=False, repr=False, compare=False)
    registry: Dict[str, Union['OperatorDict', 'UnaryOperatorDict', Callable]] = field(init=False, repr=False, compare=False, default_factory=dict) # Allow Callables for register
    numspace: Dict[str, Any] = field(init=False, repr=False, compare=False, default_factory=dict)

    # Operator dictionaries (Definitions unchanged)
    gp: 'OperatorDict' = operation_field()
    sw: 'OperatorDict' = operation_field()
    cp: 'OperatorDict' = operation_field()
    acp: 'OperatorDict' = operation_field()
    ip: 'OperatorDict' = operation_field()
    sp: 'OperatorDict' = operation_field()
    lc: 'OperatorDict' = operation_field()
    rc: 'OperatorDict' = operation_field()
    op: 'OperatorDict' = operation_field()
    rp: 'OperatorDict' = operation_field()
    proj: 'OperatorDict' = operation_field()
    add: 'OperatorDict' = operation_field()
    sub: 'OperatorDict' = operation_field()
    div: 'OperatorDict' = operation_field()
    inv: 'UnaryOperatorDict' = operation_field()
    neg: 'UnaryOperatorDict' = operation_field()
    reverse: 'UnaryOperatorDict' = operation_field()
    involute: 'UnaryOperatorDict' = operation_field()
    conjugate: 'UnaryOperatorDict' = operation_field()
    sqrt: 'UnaryOperatorDict' = operation_field()
    polarity: 'UnaryOperatorDict' = operation_field()
    unpolarity: 'UnaryOperatorDict' = operation_field()
    hodge: 'UnaryOperatorDict' = operation_field()
    unhodge: 'UnaryOperatorDict' = operation_field()
    normsq: 'UnaryOperatorDict' = operation_field()
    outerexp: 'UnaryOperatorDict' = operation_field()
    outersin: 'UnaryOperatorDict' = operation_field()
    outercos: 'UnaryOperatorDict' = operation_field()
    outertan: 'UnaryOperatorDict' = operation_field()

    # --- Metadata for Operators (Unchanged) ---
    _operator_metadata: ClassVar[Dict[str, Dict[str, Any]]] = {
        'gp': {'codegen': 'codegen_gp', 'cls': 'OperatorDict'}, 'sw': {'codegen': 'codegen_sw', 'cls': 'OperatorDict'},
        'cp': {'codegen': 'codegen_cp', 'cls': 'OperatorDict'}, 'acp': {'codegen': 'codegen_acp', 'cls': 'OperatorDict'},
        'ip': {'codegen': 'codegen_ip', 'cls': 'OperatorDict'}, 'sp': {'codegen': 'codegen_sp', 'cls': 'OperatorDict'},
        'lc': {'codegen': 'codegen_lc', 'cls': 'OperatorDict'}, 'rc': {'codegen': 'codegen_rc', 'cls': 'OperatorDict'},
        'op': {'codegen': 'codegen_op', 'cls': 'OperatorDict'}, 'rp': {'codegen': 'codegen_rp', 'cls': 'OperatorDict'},
        'proj': {'codegen': 'codegen_proj', 'cls': 'OperatorDict'}, 'add': {'codegen': 'codegen_add', 'cls': 'OperatorDict'},
        'sub': {'codegen': 'codegen_sub', 'cls': 'OperatorDict'}, 'div': {'codegen': 'codegen_div', 'cls': 'OperatorDict'},
        'inv': {'codegen': 'codegen_inv', 'cls': 'UnaryOperatorDict'}, 'neg': {'codegen': 'codegen_neg', 'cls': 'UnaryOperatorDict'},
        'reverse': {'codegen': 'codegen_reverse', 'cls': 'UnaryOperatorDict'}, 'involute': {'codegen': 'codegen_involute', 'cls': 'UnaryOperatorDict'},
        'conjugate': {'codegen': 'codegen_conjugate', 'cls': 'UnaryOperatorDict'}, 'sqrt': {'codegen': 'codegen_sqrt', 'cls': 'UnaryOperatorDict'},
        'polarity': {'codegen': 'codegen_polarity', 'cls': 'UnaryOperatorDict'}, 'unpolarity': {'codegen': 'codegen_unpolarity', 'cls': 'UnaryOperatorDict'},
        'hodge': {'codegen': 'codegen_hodge', 'cls': 'UnaryOperatorDict'}, 'unhodge': {'codegen': 'codegen_unhodge', 'cls': 'UnaryOperatorDict'},
        'normsq': {'codegen': 'codegen_normsq', 'cls': 'UnaryOperatorDict'}, 'outerexp': {'codegen': 'codegen_outerexp', 'cls': 'UnaryOperatorDict'},
        'outersin': {'codegen': 'codegen_outersin', 'cls': 'UnaryOperatorDict'}, 'outercos': {'codegen': 'codegen_outercos', 'cls': 'UnaryOperatorDict'},
        'outertan': {'codegen': 'codegen_outertan', 'cls': 'UnaryOperatorDict'},
    }

    def __post_init__(self) -> None:
        # Imports moved here to avoid potential circularity at module level
        from kingdon.multivector import MultiVector
        from kingdon.operator_dict import BladeDict

        self.numspace = {'np': np, 'MultiVector': MultiVector, 'cos': np.cos, 'sin': np.sin, 'sqrt': np.sqrt,
                         'math': math} # Added math
=======
    basis: List[str] = field(default_factory=list)
    
    # Internal attributes
    d: int = field(init=False)
    cse: bool = True
    graded: bool = False
    simp_func: Callable[[Any], Any] = field(default=lambda x: sympy.simplify(sympy.expand(x)) if isinstance(x, sympy.Expr) else x)
    wrapper: Optional[Callable] = None
    codegen_symbolcls: Any = field(default=None)
    pretty_blade: str = '𝐞'

    # Define indices for each grade
    indices_for_grades: Dict[Tuple[int, ...], List[int]] = field(init=False, default_factory=dict)

    # Operator dictionaries
    gp: OperatorDict = operation_field(metadata={'codegen': codegen_gp})
    sw: OperatorDict = operation_field(metadata={'codegen': codegen_sw})
    cp: OperatorDict = operation_field(metadata={'codegen': codegen_cp})
    acp: OperatorDict = operation_field(metadata={'codegen': codegen_acp})
    ip: OperatorDict = operation_field(metadata={'codegen': codegen_ip})
    sp: OperatorDict = operation_field(metadata={'codegen': codegen_sp})
    lc: OperatorDict = operation_field(metadata={'codegen': codegen_lc})
    rc: OperatorDict = operation_field(metadata={'codegen': codegen_rc})
    op: OperatorDict = operation_field(metadata={'codegen': codegen_op})
    rp: OperatorDict = operation_field(metadata={'codegen': codegen_rp})
    proj: OperatorDict = operation_field(metadata={'codegen': codegen_proj})
    add: OperatorDict = operation_field(metadata={'codegen': codegen_add})
    sub: OperatorDict = operation_field(metadata={'codegen': codegen_sub})
    div: OperatorDict = operation_field(metadata={'codegen': codegen_div})
    inv: UnaryOperatorDict = operation_field(metadata={'codegen': codegen_inv})
    neg: UnaryOperatorDict = operation_field(metadata={'codegen': codegen_neg})
    reverse: UnaryOperatorDict = operation_field(metadata={'codegen': codegen_reverse})
    involute: UnaryOperatorDict = operation_field(metadata={'codegen': codegen_involute})
    conjugate: UnaryOperatorDict = operation_field(metadata={'codegen': codegen_conjugate})
    sqrt: UnaryOperatorDict = operation_field(metadata={'codegen': codegen_sqrt})
    polarity: UnaryOperatorDict = operation_field(metadata={'codegen': codegen_polarity})
    unpolarity: UnaryOperatorDict = operation_field(metadata={'codegen': codegen_unpolarity})
    hodge: UnaryOperatorDict = operation_field(metadata={'codegen': codegen_hodge})
    unhodge: UnaryOperatorDict = operation_field(metadata={'codegen': codegen_unhodge})
    normsq: UnaryOperatorDict = operation_field(metadata={'codegen': codegen_normsq})
    outerexp: UnaryOperatorDict = operation_field(metadata={'codegen': codegen_outerexp})
    outersin: UnaryOperatorDict = operation_field(metadata={'codegen': codegen_outersin})
    outercos: UnaryOperatorDict = operation_field(metadata={'codegen': codegen_outercos})
    outertan: UnaryOperatorDict = operation_field(metadata={'codegen': codegen_outertan})

    def __post_init__(self) -> None:
        """Initialize the algebra with computed attributes."""
        # Import here to avoid circular imports at runtime
        from kingdon.multivector import MultiVector
        
        # Set the base namespace for compiled functions
        self.numspace = {'np': np, 'MultiVector': MultiVector}
        
        # Signature and dimension setup methods
>>>>>>> 5d1776bf2c10cd778122088d9664260375c22a81
        self._setup_signature()
        self._setup_basis_mappings()
        self._precompute_grade_indices()
        self._compute_multiplication_signs()
<<<<<<< HEAD
        self._initialize_blade_dictionary(BladeDict, MultiVector) # Sets self.pss
        self._create_registry()
        # Precompute internal gp function for cayley table
        self._gp_func = self._compile_gp_func()


    # --- Setup Methods ---
    # (No changes needed in setup methods for this step)
    def _setup_signature(self) -> None:
        if self.signature is not None:
            if not isinstance(self.signature, (list, tuple, np.ndarray)): raise TypeError("Explicit signature must be a list, tuple, or numpy array.")
            sig_list = list(self.signature); counts = Counter(sig_list)
            self.p, self.q, self.r = counts[1], counts[-1], counts[0]
            if self.p + self.q + self.r != len(sig_list): raise ValueError("Invalid signature: must contain only 1, -1, or 0")
            self.signature = np.array(sig_list)
        else:
            sig_list = [1] * self.p + [-1] * self.q + [0] * self.r; self.signature = np.array(sig_list)
        self.d = self.p + self.q + self.r
        if self.start_index is None: self.start_index = 0 if self.r == 1 else 1
        if self.basis_names is not None and len(self.basis_names) != self.d: raise ValueError(f"basis_names length ({len(self.basis_names)}) must match algebra dimension ({self.d})")

    def _setup_basis_mappings(self) -> None:
        self.canon2bin = {}; self.bin2canon = {}; self._bin2canon_prettystr = {}
        self.canon2bin["e"] = 0; self.bin2canon[0] = "e"; self._bin2canon_prettystr[0] = "1"
        vector_names = []; unicode_map = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
        for i in range(self.d):
            idx = i + self.start_index; default_name = f"e{idx}"; name = self.basis_names[i] if self.basis_names else default_name
            if name in self.canon2bin: raise ValueError(f"Duplicate basis name detected: '{name}'")
            key = 1 << i; self.canon2bin[name] = key
            if name != default_name: self.canon2bin[default_name] = key
            self.bin2canon[key] = name; vector_names.append(name)
            subscript = str(idx); pretty_sub = subscript.translate(unicode_map); self._bin2canon_prettystr[key] = f"{self.pretty_blade}{pretty_sub}"
        for grade in range(2, self.d + 1):
            for combo_indices in combinations(range(self.d), grade):
                key = sum(1 << i for i in combo_indices); sorted_canon_indices = sorted([i + self.start_index for i in combo_indices])
                default_canon_name = f"e{''.join(map(str, sorted_canon_indices))}"; self.canon2bin[default_canon_name] = key; self.bin2canon[key] = default_canon_name
                subscript = "".join(map(str, sorted_canon_indices)); pretty_sub = subscript.translate(unicode_map); self._bin2canon_prettystr[key] = f"{self.pretty_blade}{pretty_sub}"
        if self.start_index == 0: self.canon2bin["e0"] = 0

    def _blade2canon(self, basis_blade: str) -> Tuple[str, int]:
        if basis_blade == "e" or (self.start_index == 0 and basis_blade == "e0"): return "e", 0
        match_default = re.match(r'^e([0-9]+)$', basis_blade)
        if match_default:
            indices_str = match_default.group(1)
            try: indices = [int(c) for c in indices_str]
            except ValueError as e: raise ValueError(f"Invalid index character in blade string: {basis_blade}. {e}")
            if len(set(indices)) != len(indices): raise ValueError(f"Duplicate indices in blade name: {basis_blade}")
            min_idx, max_idx = self.start_index, self.start_index + self.d - 1
            if not all(min_idx <= i <= max_idx for i in indices): raise ValueError(f"Indices in '{basis_blade}' must be between {min_idx} and {max_idx}")
            sorted_indices = sorted(indices); canonical_name = f"e{''.join(map(str, sorted_indices))}"
            current_indices = list(indices); swaps = 0; n = len(current_indices)
            for i in range(n):
                for j in range(n - 1 - i):
                    if current_indices[j] > current_indices[j+1]: current_indices[j], current_indices[j+1] = current_indices[j+1], current_indices[j]; swaps += 1
            return canonical_name, swaps
        elif self.basis_names:
            matched_names = []
            remaining_blade = basis_blade
            sorted_basis_names = sorted(self.basis_names, key=len, reverse=True)
            while remaining_blade:
                found_match = False
                for name in sorted_basis_names:
                    if remaining_blade.startswith(name):
                        matched_names.append(name)
                        remaining_blade = remaining_blade[len(name):]
                        found_match = True
                        break
                if not found_match:
                    raise ValueError(f"Invalid blade string '{basis_blade}': Contains characters/substrings not matching custom basis names {self.basis_names}")
            if len(set(matched_names)) != len(matched_names): raise ValueError(f"Duplicate basis names in blade string: {basis_blade}")
            name_to_index = {name: i for i, name in enumerate(self.basis_names)}
            sorted_matched_names = sorted(matched_names, key=lambda name: name_to_index[name])
            canonical_name = "".join(sorted_matched_names)
            current_names = list(matched_names)
            swaps = 0
            n = len(current_names)
            for i in range(n):
                for j in range(n - 1 - i):
                    if name_to_index[current_names[j]] > name_to_index[current_names[j+1]]:
                        current_names[j], current_names[j+1] = current_names[j+1], current_names[j]
                        swaps += 1
            return canonical_name, swaps
        else:
             if basis_blade in self.canon2bin: return basis_blade, 0
             raise ValueError(f"Invalid basis blade format: '{basis_blade}'. Use 'e' + indices or provide custom basis_names.")

    def _precompute_grade_indices(self) -> None:
        self.indices_for_grades = defaultdict(list)
        all_keys = list(range(2**self.d))
        key_func = lambda k: bin(k).count('1')
        for grade, group in groupby(sorted(all_keys, key=key_func), key=key_func):
            self.indices_for_grades[(grade,)] = list(group)

    def _compute_multiplication_signs(self) -> None:
        self.signs = {};
        for i in range(2**self.d):
            for j in range(2**self.d):
                if i == 0 or j == 0: self.signs[(i, j)] = 1; continue
                i_indices = {k for k in range(self.d) if (i >> k) & 1}; j_indices = {k for k in range(self.d) if (j >> k) & 1}
                common_indices = i_indices.intersection(j_indices); i_only_indices = i_indices.difference(common_indices); j_only_indices = j_indices.difference(common_indices)
                swaps = 0
                for i_idx in i_only_indices:
                    for j_idx in j_only_indices:
                        if i_idx > j_idx: swaps += 1
                sign = (-1)**swaps
                for k in common_indices:
                     sig_val = self.signature[k]
                     if sig_val == 0: sign = 0; break
                     sign *= sig_val
                if sign != 0: self.signs[(i, j)] = int(sign)

    def _initialize_blade_dictionary(self, BladeDictClass: Type['BladeDict'], MVClass: Type['MultiVector']) -> None:
        lazy_init = self.d > 6; self.blades = BladeDictClass(self, lazy=lazy_init)
        if self.d > 0:
            pseudoscalar_key = (1 << self.d) - 1; pseudoscalar_name = self.bin2canon.get(pseudoscalar_key)
            if pseudoscalar_name: self.pss = self.blades[pseudoscalar_name]
            else:
                log.warning(f"Could not find canonical name for pseudoscalar key {pseudoscalar_key}. Creating manually.")
                self.pss = MVClass.fromkeysvalues(self, keys=(pseudoscalar_key,), values=[1])
        else:
            self.pss = MVClass.fromkeysvalues(self, keys=(0,), values=[1])

    def _create_registry(self) -> None:
        from kingdon.operator_dict import OperatorDict, UnaryOperatorDict
        import kingdon.codegen as codegen_module
        self.registry = {}
        for op_name, meta in self._operator_metadata.items():
            codegen_func_name = meta.get('codegen'); operator_cls_name = meta.get('cls', 'OperatorDict')
            if codegen_func_name:
                codegen_func = getattr(codegen_module, codegen_func_name, None)
                if codegen_func is None: log.warning(f"Codegen function '{codegen_func_name}' not found for operator '{op_name}'. Skipping."); continue
                operator_cls = UnaryOperatorDict if operator_cls_name == 'UnaryOperatorDict' else OperatorDict
                operator_instance = operator_cls(name=op_name, algebra=self, codegen=codegen_func)
                self.registry[op_name] = operator_instance; setattr(self, op_name, operator_instance)

    def _compile_gp_func(self) -> Callable[[int, int], Tuple[int, int]]:
        @lru_cache(maxsize=None)
        def gp_basis(b1_idx: int, b2_idx: int) -> Tuple[int, int]:
            """ Computes basis_out_idx, sign = gp(basis[b1_idx], basis[b2_idx])."""
            sign = self.signs.get((b1_idx, b2_idx), 0)
            result_idx = b1_idx ^ b2_idx
            return result_idx, sign
        return gp_basis

    # ========================================================================
    # Helper Methods
    # ========================================================================

    def _parse_composite_custom_key(self, key_str: str) -> Optional[List[int]]:
        # (Unchanged)
        if not self.basis_names:
            return None # Cannot parse if no custom names defined
        matched_keys = []
        matched_names_set = set()
        remaining_key_str = key_str
        sorted_basis_names = sorted(self.basis_names, key=len, reverse=True)
        while remaining_key_str:
            found_match = False
            for name in sorted_basis_names:
                if remaining_key_str.startswith(name):
                    int_key = self.canon2bin.get(name)
                    if int_key is None:
                        log.error(f"Internal error: Custom basis name '{name}' not found in canon2bin.")
                        return None
                    if name in matched_names_set:
                        log.debug(f"Parsing failed: Duplicate component '{name}' found in '{key_str}'.")
                        return None
                    matched_keys.append(int_key)
                    matched_names_set.add(name)
                    remaining_key_str = remaining_key_str[len(name):]
                    found_match = True
                    break
            if not found_match:
                log.debug(f"Parsing failed: Cannot match remaining part '{remaining_key_str}' of '{key_str}' with custom basis names.")
                return None
        return matched_keys


    # ========================================================================
    # Multivector Creation Methods (Using _process_creation_args)
    # ========================================================================

    def _process_creation_args(self, *args: Any, **kwargs: Any) -> Tuple[Optional[Tuple[int, ...]], Optional[Union[List[Any], np.ndarray, Dict[str, Any]]], Optional[Tuple[int, ...]], Optional[str], Optional[Any]]:
        """Internal helper to parse arguments for multivector creation."""
        values = kwargs.pop('values', None); keys = kwargs.pop('keys', None); grades = kwargs.pop('grades', None); name = kwargs.pop('name', None); symbolcls = kwargs.pop('symbolcls', None)
        # --- Argument Conflict Checks ---
        if keys is not None and grades is not None:
            raise ValueError("Cannot specify both 'keys' and 'grades' arguments.")
        if keys is not None and kwargs:
            raise ValueError("Cannot specify both 'keys' and keyword arguments for blades.")
        if grades is not None and kwargs:
             raise ValueError("Cannot specify 'grades' and keyword arguments for blades. Use 'values' list/array with 'grades'.")
        if name is not None and kwargs:
             raise ValueError("Cannot specify both 'name' (symbolic) and keyword arguments for blades.")
        if name is not None and values is not None:
             raise ValueError("Cannot specify both 'name' (symbolic) and explicit 'values'.")

        # --- Process Positional Args ---
        if args:
            if len(args) > 1: raise TypeError(f"Too many positional arguments ({len(args)}), expected 0 or 1.")
            arg = args[0]
            if isinstance(arg, str) and name is None: name = arg
            elif values is None: values = arg
            # Allow positional scalar value if no name/values set yet
            elif isinstance(arg, (int, float, complex, sympy.Expr)) and values is None and name is None:
                values = arg
            else:
                raise TypeError(f"Unexpected positional argument type '{type(arg).__name__}' or argument already specified via keyword.")

        # --- Process Keyword Args (Blades) ---
        if kwargs:
             if not hasattr(self, 'canon2bin') or not self.canon2bin: self._setup_basis_mappings()
             blade_kwargs = {}
             for k, v in kwargs.items():
                  # Check if k is a valid blade name (standard, custom, or composite custom)
                  is_valid_blade_key = False
                  if isinstance(k, str):
                       try:
                           self._blade2canon(k) # Test if it's a standard or custom name
                           is_valid_blade_key = True
                       except ValueError:
                           # Try parsing as composite custom key
                           if self._parse_composite_custom_key(k):
                               is_valid_blade_key = True
                  elif k == 'e': # Allow 'e' specifically
                       is_valid_blade_key = True

                  if is_valid_blade_key:
                      blade_kwargs[k] = v
                  else:
                      raise TypeError(f"{type(self).__name__}.multivector() got an unexpected keyword argument '{k}'. Valid blade names/keys are expected.")

             if blade_kwargs:
                 if values is not None and not isinstance(values, dict):
                      raise TypeError("Cannot specify both 'values' (as list/array) and keyword arguments for blades.")
                 if values is None: values = {}
                 elif not isinstance(values, dict):
                      raise TypeError("'values' must be a dictionary when using keyword arguments for blades.")
                 # Check for overlaps between dict values and kwargs
                 overlapping_keys = set(values.keys()) & set(blade_kwargs.keys())
                 if overlapping_keys:
                      raise ValueError(f"Blade keys specified in both 'values' dict and keyword arguments: {overlapping_keys}")
                 values.update(blade_kwargs)

        return keys, values, grades, name, symbolcls


    def multivector(self, *args: Any, **kwargs: Any) -> 'MultiVector':
        """ Primary factory method to create a multivector in this algebra. """
        # (Docstring from previous step retained)
        from kingdon.multivector import MultiVector
        keys, values, grades, name, symbolcls_arg = self._process_creation_args(*args, **kwargs)
        symbolcls_arg = symbolcls_arg or self.codegen_symbolcls

        # --- Variable Initialization ---
        final_keys: Optional[Tuple[int, ...]] = None
        final_values: Union[List[Any], np.ndarray]
        target_len: int = -1 # Expected length based on keys/grades

        # --- Determine Final Keys and Target Length ---
        max_key_val = (1 << self.d) - 1 # Maximum valid integer key

        if keys is not None: # Case 1: Explicit keys
            processed_keys = []
            keys_input = (keys,) if isinstance(keys, (int, str)) else keys # Normalize
            if not isinstance(keys_input, (tuple, list)):
                raise TypeError(f"`keys` argument must be an int, string, or a tuple/list. Got {type(keys_input).__name__}.")

            for k in keys_input:
                if isinstance(k, int):
                    if not (0 <= k <= max_key_val):
                         raise ValueError(f"Integer key {k} out of range [0, {max_key_val}] for dimension {self.d}.")
                    processed_keys.append(k)
                elif isinstance(k, str):
                    try:
                        canon_name, _ = self._blade2canon(k) # Validates format/indices
                        int_key = self.canon2bin.get(canon_name)
                        if int_key is None: raise KeyError(f"Internal error mapping canonical name '{canon_name}'.")
                        processed_keys.append(int_key)
                    except (ValueError, KeyError) as e:
                        raise ValueError(f"Invalid basis blade name '{k}' provided in 'keys': {e}") from e
                else:
                    raise TypeError(f"Items in `keys` must be integers or strings. Got type {type(k).__name__}.")
            final_keys = tuple(sorted(list(set(processed_keys)), key=self._canonical_key_sort_func())) # Use set for uniqueness
            target_len = len(final_keys)

        elif grades is not None: # Case 2: Grades provided
            grades_input = (grades,) if isinstance(grades, int) else grades # Normalize
            if not isinstance(grades_input, (tuple, list)):
                raise TypeError(f"`grades` must be an int or a tuple/list. Got {type(grades_input).__name__}.")
            valid_grades = set()
            for g in grades_input:
                 if not isinstance(g, int) or not (0 <= g <= self.d):
                      raise ValueError(f"Grade {g} is invalid. Grades must be integers between 0 and {self.d}.")
                 valid_grades.add(g)

            collected_keys = []
            for g in sorted(list(valid_grades)):
                grade_keys = self.indices_for_grades.get((g,), [])
                collected_keys.extend(grade_keys)
            final_keys = tuple(sorted(collected_keys, key=self._canonical_key_sort_func()))
            target_len = len(final_keys)
            if not final_keys and self.d == 0 and 0 in valid_grades: final_keys = (0,); target_len = 1

        elif isinstance(values, dict): # Case 3: Values as Dict/Keywords
            processed_keys_dict = {}
            for k_in, v in values.items(): # Use values dict directly (already includes kwargs)
                 int_key: Optional[int] = None
                 if isinstance(k_in, str):
                     try:
                         # Try standard interpretation first
                         canon_name, _ = self._blade2canon(k_in)
                         int_key_opt = self.canon2bin.get(canon_name)
                         # If standard name fails, try parsing as composite custom key
                         if int_key_opt is None:
                              parsed_custom_keys = self._parse_composite_custom_key(k_in)
                              if parsed_custom_keys is not None:
                                   combined_key = reduce(operator.xor, parsed_custom_keys, 0)
                                   int_key_opt = combined_key
                              else: raise KeyError # Parsing failed
                         int_key = int_key_opt
                     except (ValueError, KeyError) as e:
                         raise KeyError(f"Invalid basis blade key '{k_in}': {e}") from e
                 elif isinstance(k_in, int):
                     if not (0 <= k_in <= max_key_val):
                          raise KeyError(f"Integer key {k_in} out of range [0, {max_key_val}].")
                     int_key = k_in
                 else:
                     raise TypeError(f"Invalid key type in values dictionary: {type(k_in).__name__}. Use int or str.")

                 if int_key is None: raise RuntimeError(f"Internal error: Failed to determine integer key for input key '{k_in}'.")
                 if int_key in processed_keys_dict:
                     raise ValueError(f"Duplicate key specified: '{k_in}' resolves to the same component as '{processed_keys_dict[int_key]}'.")
                 processed_keys_dict[int_key] = k_in

            final_keys = tuple(sorted(processed_keys_dict.keys(), key=self._canonical_key_sort_func()))
            target_len = len(final_keys)
            if not final_keys: final_keys = (0,); target_len = 1;

        elif name is not None: # Case 4: Symbolic
             # Determine keys based on whether grades were also specified implicitly
             if grades is not None: # Symbolic + Grades (handled in Case 2)
                 pass # final_keys/target_len already set
             else: # Full symbolic MV
                 final_keys = tuple(k for g_keys in sorted(self.indices_for_grades.values()) for k in g_keys)
                 target_len = len(self) # Length of full algebra
                 if not final_keys and self.d == 0: final_keys = (0,); target_len = 1

        else: # Case 5: Values as list/array/scalar or None (default)
            if values is None: final_keys, target_len = (0,), 1
            elif isinstance(values, (int, float, complex, sympy.Expr)): final_keys, target_len = (0,), 1
            elif hasattr(values, '__len__'):
                val_len = len(values)
                if val_len == len(self): # Assume full canonical sequence
                    final_keys = tuple(k for g_keys in sorted(self.indices_for_grades.values()) for k in g_keys)
                    target_len = len(self)
                    if not final_keys and self.d == 0: final_keys = (0,); target_len = 1
                else: # Ambiguous length - ERROR
                     raise ValueError(f"Length of sequence 'values' ({val_len}) is ambiguous. Provide explicit 'keys' or 'grades', or ensure length matches algebra size ({len(self)}).")
            else: final_keys, target_len = (0,), 1 # Default for unknown type

        # --- Final Key/Length Consistency Check ---
        if final_keys is None or target_len == -1:
             raise RuntimeError("Internal state error: final_keys or target_len could not be determined.")

        # --- Determine Final Values ---
        default_zero = sympy.Integer(0) if self.codegen_symbolcls != float else 0.0

        if name is not None: # Case 4: Symbolic
             sym_func = symbolcls_arg or sympy.Symbol; final_values_dict = {}
             for k in final_keys:
                 blade_name_suffix = self.bin2canon.get(k, f'k{k}')
                 sym_name = f"{name}_{blade_name_suffix}"
                 final_values_dict[k] = default_zero if k == 0 else sym_func(sym_name) # Use appropriate zero
             final_values = [final_values_dict[k] for k in final_keys]

        elif isinstance(values, dict): # Case 3: Dict/Keywords
             original_key_to_int = {v: k for k, v in processed_keys_dict.items()}
             final_values_map = {}
             for k_orig, v in values.items():
                 int_key = original_key_to_int.get(k_orig)
                 if int_key is None: raise RuntimeError(f"Internal error mapping original key '{k_orig}'")
                 final_values_map[int_key] = sympify(v) if isinstance(v, (int, float, complex)) and self.codegen_symbolcls != float else v
             final_values = [final_values_map.get(k, default_zero) for k in final_keys]
             if not values and final_keys == (0,): final_values = [default_zero]

        elif isinstance(values, (list, tuple, np.ndarray)): # Case 1, 2, 5 (Sequence)
             is_empty = (isinstance(values, np.ndarray) and values.size == 0) or (not isinstance(values, np.ndarray) and not values)
             if is_empty:
                 if target_len > 0: final_values = [default_zero] * target_len
                 else: final_values = []
             elif len(values) != target_len:
                  # Improved error message for sequence length mismatch
                  spec_type = "explicit keys" if keys is not None else ("specified grades" if grades is not None else "full algebra")
                  raise ValueError(f"Length of sequence 'values' ({len(values)}) does not match the expected length ({target_len}) for the {spec_type}.")
             else: final_values = list(values) if not isinstance(values, np.ndarray) else values

        elif isinstance(values, (int, float, complex, sympy.Expr)): # Case 5 (Scalar)
             if final_keys == (0,) and target_len == 1: final_values = [values]
             else: raise ValueError(f"Scalar value '{values}' provided, but target keys {final_keys} indicate a non-scalar structure.")

        elif values is None: # Case 5 (Default) or Case 2 (Grades only)
             final_values = [default_zero] * target_len
             if final_keys == (0,) and target_len == 1: final_values = [default_zero]

        else: raise TypeError(f"Unsupported type for 'values': {type(values).__name__}")

        # --- Final Construction ---
        final_keys_tuple = tuple(final_keys) if not isinstance(final_keys, tuple) else final_keys
        return MultiVector.fromkeysvalues(self, keys=final_keys_tuple, values=final_values)

    # ========================================================================
    # Convenience Factory Methods (Calling self.multivector)
    # ========================================================================
    # (Methods scalar, vector, bivector, trivector, purevector, evenmv, oddmv
    #  blade, pseudoscalar, rotor, translator, reflector, motor remain unchanged
    #  from previous step, as they already call self.multivector)
    def scalar(self, value: Any = 0, *, name: Optional[str] = None) -> 'MultiVector':
        return self.multivector(values=[value], grades=(0,), name=name)

    def vector(self, values: Union[Sequence[Any], np.ndarray, Dict[Union[str, int], Any]], *, name: Optional[str] = None) -> 'MultiVector':
        return self.multivector(values=values, grades=(1,), name=name)

    def bivector(self, values: Union[Sequence[Any], np.ndarray, Dict[Union[str, int], Any]], *, name: Optional[str] = None) -> 'MultiVector':
        return self.multivector(values=values, grades=(2,), name=name)

    def trivector(self, values: Union[Sequence[Any], np.ndarray, Dict[Union[str, int], Any]], *, name: Optional[str] = None) -> 'MultiVector':
        return self.multivector(values=values, grades=(3,), name=name)

    def purevector(self, values: Union[Sequence[Any], np.ndarray, Dict[Union[str, int], Any]], grade: int, *, name: Optional[str] = None) -> 'MultiVector':
        if not isinstance(grade, int) or not (0 <= grade <= self.d):
             raise ValueError(f"Grade must be an integer between 0 and {self.d}")
        return self.multivector(values=values, grades=(grade,), name=name)

    def evenmv(self, values: Union[Sequence[Any], np.ndarray, Dict[Union[str, int], Any]], *, name: Optional[str] = None) -> 'MultiVector':
        even_grades = tuple(g for g in range(0, self.d + 1) if g % 2 == 0)
        if not even_grades and self.d == 0: even_grades = (0,)
        elif not even_grades: return self.multivector(values={}, grades=())
        return self.multivector(values=values, grades=even_grades, name=name)

    def oddmv(self, values: Union[Sequence[Any], np.ndarray, Dict[Union[str, int], Any]], *, name: Optional[str] = None) -> 'MultiVector':
        odd_grades = tuple(g for g in range(1, self.d + 1) if g % 2 != 0)
        if not odd_grades: return self.multivector(values={}, grades=())
        return self.multivector(values=values, grades=odd_grades, name=name)

    def blade(self, spec: Union[str, Sequence[int]], value: Any = 1, *, name: Optional[str] = None) -> 'MultiVector':
        if isinstance(spec, (list, tuple)):
            if not spec: return self.scalar(value, name=name)
            min_idx, max_idx = self.start_index, self.start_index + self.d - 1
            if not all(isinstance(i, int) and min_idx <= i <= max_idx for i in spec):
                 raise ValueError(f"Indices {spec} must be integers between {min_idx} and {max_idx}")
            if len(set(spec)) != len(spec): raise ValueError(f"Indices in {spec} must be unique.")
            sorted_indices = sorted(list(spec))
            key_str = f"e{''.join(map(str, sorted_indices))}"
            current_indices = list(spec); swaps = 0; n = len(current_indices)
            for i in range(n):
                 for j in range(n - 1 - i):
                      if current_indices[j] > current_indices[j+1]: current_indices[j], current_indices[j+1] = current_indices[j+1], current_indices[j]; swaps += 1
            sign = (-1)**swaps; final_value = sign * value
        elif isinstance(spec, str):
            try:
                canonical_name, swaps = self._blade2canon(spec)
                key_str = canonical_name
                sign = (-1)**swaps; final_value = sign * value
            except ValueError as e: raise ValueError(f"Invalid blade specification string '{spec}': {e}")
        else: raise TypeError("Blade specification must be a string or a list/tuple of indices.")
        try: return self.multivector(values={key_str: final_value}, name=name)
        except KeyError as e: raise ValueError(f"Could not create blade for spec '{spec}': {e}")

    def pseudoscalar(self, value: Any = 1, *, name: Optional[str] = None) -> 'MultiVector':
        if value == 1 and name is None: return self.pss
        else: scaled_pss = self.pss * value; #if name: scaled_pss.name = name; # Name assignment problematic on numeric MV
        return scaled_pss

    def rotor(self, angle: float, plane_indices: Tuple[int, int], *, name: Optional[str] = None) -> 'MultiVector':
        if not isinstance(plane_indices, tuple) or len(plane_indices) != 2: raise TypeError("plane_indices must be a tuple of two integers.")
        i, j = sorted(plane_indices)
        if i == j: raise ValueError("plane_indices must contain two distinct indices.")
        scalar_part = self.scalar(math.cos(angle / 2.0)); bivec_part = self.blade([i, j], value=math.sin(angle / 2.0))
        rotor = scalar_part + bivec_part; #if name: rotor.name = name;
        return rotor

    def translator(self, direction: Sequence[float], distance: float = 1.0, *, name: Optional[str] = None) -> 'MultiVector':
        if self.r != 1 or self.start_index != 0: raise ValueError("Translator requires a PGA setup (r=1, start_index=0).")
        if len(direction) != self.d - 1: raise ValueError(f"Direction vector must have {self.d - 1} components for {self.d}D PGA.")
        translator = self.scalar(1)
        for i, component in enumerate(direction, start=1):
            if component != 0:
                coeff = (distance / 2.0) * component; blade_term = self.blade([0, i], value=coeff)
                translator += blade_term
        #if name: translator.name = name;
        return translator

    def reflector(self, normal: Sequence[float], *, name: Optional[str] = None) -> 'MultiVector':
        if len(normal) != self.d: raise ValueError(f"Normal vector length ({len(normal)}) must match algebra dimension ({self.d}).")
        norm_vec = self.vector(values=normal, name=name)
        try: reflector = norm_vec.normalized(); #if name: reflector.name = name;
        except ZeroDivisionError: raise ValueError("Normal vector cannot be zero (cannot normalize).")
        except TypeError as e: raise TypeError(f"Could not normalize the normal vector: {e}")
        return reflector

    def motor(self, rotation: Optional[Tuple[float, Tuple[int, int]]] = None, translation: Optional[Tuple[Sequence[float], float]] = None, *, name: Optional[str] = None) -> 'MultiVector':
        if not rotation and not translation: return self.scalar(1, name=name)
        if translation and (self.r != 1 or self.start_index != 0): raise ValueError("Translation requires a PGA setup (r=1, start_index=0).")
        R = self.scalar(1); T = self.scalar(1)
        if rotation: angle, plane_indices = rotation; R = self.rotor(angle, plane_indices)
        if translation: direction, distance = translation; T = self.translator(direction, distance)
        motor = T * R; #if name: motor.name = name;
        return motor


    # ========================================================================
    # Properties and Other Methods (Existing)
    # ========================================================================
    # (Methods pss_inv, cayley, frame, reciprocal_frame, register, graph,
    # __len__, __eq__, __hash__, grade, dual, undual, from_signature, from_name
    # remain unchanged from previous step)
    @property
    @lru_cache(maxsize=None)
    def pss_inv(self) -> 'MultiVector':
        if not hasattr(self, 'pss'): raise AttributeError("Algebra not fully initialized, pss attribute missing.")
        try: return self.pss.inv()
        except ZeroDivisionError as e: raise ZeroDivisionError(f"Cannot compute pss_inv: Pseudoscalar is not invertible. Original error: {e}")

    @property
    def cayley(self) -> np.ndarray:
        basis_indices = np.arange(len(self)); vectorized_gp = np.vectorize(self._gp_func, otypes=[object])
        i_indices, j_indices = np.meshgrid(basis_indices, basis_indices, indexing='ij'); table = vectorized_gp(i_indices, j_indices)
        return table

    @property
    def frame(self) -> List['MultiVector']:
        vector_keys = self.indices_for_grades.get((1,), []); return [self.blades[self.bin2canon[key]] for key in vector_keys]

    @property
    def reciprocal_frame(self) -> List['MultiVector']:
        if self.r > 0: raise NotImplementedError("Reciprocal frame is not uniquely defined for degenerate algebras (r > 0).")
        if self.d == 0: return []
        n = self.d; basis_vectors = self.frame; pss_inv_val = self.pss_inv; reciprocal_vectors = []
        all_vector_indices = list(range(n))
        for i in range(n):
            wedge_indices = [j for j in all_vector_indices if j != i]
            if not wedge_indices: wedge_product = self.scalar(1)
            else: wedge_product = reduce(operator.xor, [basis_vectors[j] for j in wedge_indices])
            sign = (-1)**i; reciprocal_vec = sign * wedge_product * pss_inv_val
            reciprocal_vectors.append(reciprocal_vec)
        return reciprocal_vectors

    def register(self, _func: Optional[Callable] = None, *, symbolic: bool = False, **options) -> Callable:
        def decorator_register(func: Callable) -> Callable:
            @wraps(func)
            def wrapper(*args, **kwargs):
                log.info(f"Calling registered function '{func.__name__}' (symbolic={symbolic}, options={options})")
                return func(*args, **kwargs)
            log.info(f"Registering function '{func.__name__}' with algebra {self}. Options: symbolic={symbolic}, {options}")
            self.registry[func.__name__] = wrapper
            return wrapper
        if _func is None: return decorator_register
        else: return decorator_register(_func)

    def graph(self, *args, **kwargs) -> Any:
        try:
            from kingdon.graph import GraphWidget
            from kingdon.multivector import MultiVector
        except ImportError: log.error("GraphWidget or MultiVector could not be imported. Please ensure 'kingdon[graph]' extras are installed."); raise
        mv_args = [arg for arg in args if isinstance(arg, MultiVector)]
        return GraphWidget(self, *mv_args, **kwargs)

    def __len__(self) -> int: return 2**self.d

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, Algebra): return NotImplemented
        sig_equal = np.array_equal(self.signature, other.signature) if self.signature is not None and other.signature is not None else self.signature == other.signature
        basis_names_equal = self.basis_names == other.basis_names
        return (self.p == other.p and self.q == other.q and self.r == other.r and self.start_index == other.start_index and sig_equal and basis_names_equal and self.graded == other.graded)

    def __hash__(self) -> int:
        sig_tuple = tuple(self.signature.tolist()) if self.signature is not None else None
        basis_tuple = tuple(self.basis_names) if self.basis_names else None
        return hash((self.p, self.q, self.r, sig_tuple, self.start_index, basis_tuple, self.graded))

    def grade(self, mv: 'MultiVector', *grades: int) -> 'MultiVector': return mv.grade(*grades)
    def dual(self, mv: 'MultiVector', kind: str = 'auto') -> 'MultiVector': return mv.dual(kind=kind)
    def undual(self, mv: 'MultiVector', kind: str = 'auto') -> 'MultiVector': return mv.undual(kind=kind)

    @classmethod
    def from_signature(cls, signature: List[int], start_index: Optional[int] = None) -> 'Algebra': return cls(signature=signature, start_index=start_index)
    @classmethod
    def from_name(cls, name: str) -> 'Algebra':
        name_upper = name.upper().strip()
        algebras = {"CGA2D": {'p': 3, 'q': 1, 'r': 0}, "CGA3D": {'p': 4, 'q': 1, 'r': 0}, "PGA2D": {'p': 2, 'q': 0, 'r': 1, 'start_index': 0}, "PGA3D": {'p': 3, 'q': 0, 'r': 1, 'start_index': 0}, "3DPGA": {'p': 3, 'q': 0, 'r': 1, 'start_index': 0}, "VGA2D": {'p': 2, 'q': 0, 'r': 0}, "VGA3D": {'p': 3, 'q': 0, 'r': 0}, "STA":   {'p': 1, 'q': 3, 'r': 0}, "APS":   {'p': 3, 'q': 0, 'r': 0}}
        if name_upper in algebras: return cls(**algebras[name_upper])
        else:
             match_vga = re.match(r'^VGA\((\d+)(?:,(\d+))?(?:,(\d+))?\)$', name_upper)
             if match_vga: p = int(match_vga.group(1)); q = int(match_vga.group(2) or 0); r = int(match_vga.group(3) or 0); return cls(p=p, q=q, r=r)
             match_cl = re.match(r'^CL\((\d+)(?:,(\d+))?(?:,(\d+))?\)$', name_upper)
             if match_cl: p = int(match_cl.group(1)); q = int(match_cl.group(2) or 0); r = int(match_cl.group(3) or 0); return cls(p=p, q=q, r=r)
             raise ValueError(f"Unknown or invalid algebra name: {name}")

# Convenience factory functions calling Algebra.from_name (Unchanged)
def PGA2D() -> Algebra: return Algebra.from_name("PGA2D")
def PGA3D() -> Algebra: return Algebra.from_name("PGA3D")
def CGA2D() -> Algebra: return Algebra.from_name("CGA2D")
def CGA3D() -> Algebra: return Algebra.from_name("CGA3D")
def VGA3D() -> Algebra: return Algebra.from_name("VGA3D")
def STA() -> Algebra: return Algebra.from_name("STA")
APS = VGA3D
=======
        self._initialize_blade_dictionary()
        self._create_registry()

    def _setup_signature(self) -> None:
        """Set up the algebra signature and dimensions."""
        if self.signature is not None:
            counts = Counter(self.signature)
            self.p, self.q, self.r = counts[1], counts[-1], counts[0]
            if self.p + self.q + self.r != len(self.signature):
                raise TypeError("Invalid signature: must contain only 1, -1, or 0")
            self.signature = np.array(self.signature)
        else:
            self.signature = np.array([1] * self.p + [-1] * self.q + [0] * self.r)
        
        self.d = self.p + self.q + self.r
        self.start_index = 0 if (self.r == 1 and self.start_index is None) else self.start_index or 1

    def _setup_basis_mappings(self) -> None:
        """Set up mappings between basis blade names and binary keys."""
        # Create the necessary attributes
        self.canon2bin = {}
        self.bin2canon = {}
        self._bin2canon_prettystr = {}
        
        # Map the scalar (identity) basis blade
        self.canon2bin["e"] = 0
        self.bin2canon[0] = "e"
        self._bin2canon_prettystr[0] = "1"
        
        # Map the remaining basis blades
        for i in range(self.d):
            # Binary representation uses 1-shifted indices
            idx = i + self.start_index
            
            # Canonical name
            name = f"e{idx}"
            self.canon2bin[name] = 1 << i
            self.bin2canon[1 << i] = name
            self._bin2canon_prettystr[1 << i] = f"{self.pretty_blade}₁" if i == 0 else f"{self.pretty_blade}_{idx}"
        
        # Generate composite basis blades
        for i in range(2, 2**self.d):
            # Skip if not a valid blade
            if i not in self.bin2canon:
                # Binary decomposition to generate name
                indices = [j + self.start_index for j in range(self.d) if (i & (1 << j))]
                
                # Create canonical name and mapping
                name = f"e{''.join(str(idx) for idx in indices)}"
                self.canon2bin[name] = i
                self.bin2canon[i] = name
                
                # Create pretty string version with subscripts
                indices_str = ''.join(str(idx) for idx in indices)
                self._bin2canon_prettystr[i] = f"{self.pretty_blade}_{{{indices_str}}}"

    def _blade2canon(self, basis_blade: str) -> Tuple[str, int]:
        """
        Convert a blade string to canonical form and determine sign.
        
        Args:
            basis_blade: The basis blade string (e.g., 'e12', 'e21')
        
        Returns:
            Tuple containing the canonical form and number of swaps needed
        
        Example:
            >>> self._blade2canon('e21')
            ('e12', 1)  # One swap needed, so sign is -1
        """
        # Handle special cases
        if basis_blade in ("e", "e0"):
            return basis_blade, 0
            
        # Extract indices from the blade string
        match = re.match(r'^e([0-9a-fA-F]+)$', basis_blade)
        if not match:
            raise ValueError(f"Invalid basis blade format: {basis_blade}")
            
        # Get the indices and sort them
        indices = [int(c) for c in match.group(1)]
        sorted_indices = sorted(indices)
        
        # Count the number of swaps needed to sort
        swaps = 0
        # Simple bubble sort to count swaps
        for i in range(len(indices)):
            for j in range(i + 1, len(indices)):
                if indices[i] > indices[j]:
                    indices[i], indices[j] = indices[j], indices[i]
                    swaps += 1
                    
        # Create the canonical form
        canonical = f"e{''.join(str(idx) for idx in sorted_indices)}"
        
        return canonical, swaps

    def _precompute_grade_indices(self) -> None:
        """Precompute indices for each grade combination."""
        # Create a mapping from grade tuple to list of indices
        self.indices_for_grades = {}
        
        # For each grade or combination of grades, collect the relevant indices
        for i in range(2**self.d):
            grade = bin(i).count('1')
            
            # Add to single grade
            if (grade,) not in self.indices_for_grades:
                self.indices_for_grades[(grade,)] = []
            self.indices_for_grades[(grade,)].append(i)
            
            # Add to combined grades (0 to grade)
            grade_tuple = tuple(range(grade + 1))
            if grade_tuple not in self.indices_for_grades:
                self.indices_for_grades[grade_tuple] = []
            self.indices_for_grades[grade_tuple].append(i)

    def _compute_multiplication_signs(self) -> None:
        """Compute multiplication signs for basis blade products."""
        # Initialize sign table
        self.signs = {}
    
        # Precompute signs for basis vector multiplication
        signs = np.zeros((self.d, self.d), dtype=int)
        for i in range(self.d):
            signs[i, i] = int(self.signature[i])
    
        # Compute signs for all basis blade combinations
        for i in range(2**self.d):
            for j in range(2**self.d):
                sign = self._compute_product_sign(i, j, signs)
                if sign != 0:
                    self.signs[(i, j)] = sign
    
        # Add additional combinations for completeness
        # This includes all valid keys that might be used in operations
        for i in range(2**self.d):
            # Scalar multiplication always has sign 1
            self.signs[(0, i)] = 1
            self.signs[(i, 0)] = 1
        
            # Self-multiplication gives signature
            if i > 0 and (i & (i - 1)) == 0:  # Check if i is a power of 2 (single basis vector)
                bit_pos = i.bit_length() - 1
                if bit_pos < len(self.signature):
                    self.signs[(i, i)] = int(self.signature[bit_pos])
                else:
                    self.signs[(i, i)] = 1  # Default for out-of-range
                    
    def _compute_product_sign(self, a: int, b: int, signs: np.ndarray) -> int:
        """Compute the sign for the product of two basis blades."""
        # Base cases
        if a == 0:
            return 1
        if b == 0:
            return 1
        
        # Decompose into basis vectors
        a_bits = [i for i in range(self.d) if (a & (1 << i))]
        b_bits = [i for i in range(self.d) if (b & (1 << i))]
        
        # Count sign changes
        sign = 1
        for i in a_bits:
            for j in b_bits:
                if i == j:
                    sign *= signs[i, j]
                elif i > j:
                    sign *= -1
        
        return sign
    
    def _initialize_blade_dictionary(self) -> None:
        self.blades = BladeDict(self, lazy=True)
        pseudoscalar_name = f"e{''.join(str(i + self.start_index) for i in range(self.d))}"
        try:
            self.pss = self.blades[pseudoscalar_name]
        except Exception as e:
            warnings.warn(f"Could not create pseudoscalar {pseudoscalar_name}: {e}", stacklevel=2)
            highest_key = (1 << self.d) - 1
            self.pss = self.multivector(values=[1], keys=(highest_key,))


    

    def _create_registry(self) -> None:
        """Create the registry of operators."""
        self.registry = {}
        for f in self.__dataclass_fields__.values():
            if f.metadata.get('codegen'):
                operator_type = (
                    f.type if isinstance(f.type, type) and 
                    issubclass(f.type, (OperatorDict, UnaryOperatorDict)) 
                    else OperatorDict
                )

                operator_instance = operator_type(
                    name=f.name, 
                    algebra=self, 
                    codegen=f.metadata['codegen']
                )
                self.registry[f.name] = operator_instance

        # Set attributes for registry items
        for name, op in self.registry.items():
            setattr(self, name, op)

    # ========================================================================
    # Multivector Creation Methods
    # ========================================================================
    
    def multivector(self, *args: Any, **kwargs: Any) -> 'MultiVector':
        """
        Create a new multivector in this algebra.
        
        Args:
            *args: Positional arguments
            **kwargs: Keyword arguments including:
                - values: Coefficients for basis blades
                - keys: Binary keys for basis blades
                - name: Optional name for symbolic representation
                - grades: Filter to specific grades
                
        Returns:
            A new multivector
            
        Example:
            >>> # Create a general multivector
            >>> mv = algebra.multivector(values=[1, 2, 3], grades=(0, 1))
            >>> # Create a named symbolic multivector
            >>> x = algebra.multivector(name="x", grades=(1,))
        """
        # Import here to avoid circular imports at runtime
        from kingdon.multivector import MultiVector
        
        # Directly create a MultiVector instance
        mv = object.__new__(MultiVector)
        mv.algebra = self
        
        # Initialize the multivector with the provided arguments
        MultiVector.__init__(mv, self, **kwargs)
        
        # Create keys and values based on the arguments
        if 'values' in kwargs and 'grades' in kwargs:
            # Get the dimension
            d = self.d
            grades = kwargs['grades']
            values = kwargs['values']
            
            # Convert to a list if it's a single value
            if not isinstance(values, (list, tuple, np.ndarray)):
                values = [values]
            
            # Determine the keys based on the grades
            if len(grades) == 1:
                grade = grades[0]
                if grade == 1:  # Vector case
                    keys = tuple(1 << i for i in range(d))
                elif grade == 2 and d >= 3 and len(values) == 3:  # Bivector case
                    keys = (3, 5, 6)
                elif grade == 0:  # Scalar case
                    keys = (0,)
                else:
                    # General case: all basis blades of the specified grade
                    keys = tuple(i for i in range(2**d) if bin(i).count('1') == grade)
            else:
                # Multiple grades: collect all basis blades of the specified grades
                keys = tuple(i for i in range(2**d) if bin(i).count('1') in grades)
            
            # Set the values and keys
            mv._keys = keys
            mv._values = values
            mv.shape = values.shape if isinstance(values, np.ndarray) and values.ndim > 1 else ()
        
        return mv

    def scalar(self, value: Union[int, float, List, Tuple, np.ndarray] = 1, 
               name: Optional[str] = None) -> 'MultiVector':
        """
        Create a scalar multivector (grade 0).
        
        Args:
            value: The scalar value (default: 1)
            name: Optional name for symbolic representation
            
        Returns:
            A new scalar multivector
            
        Example:
            >>> # Create a scalar with value 5
            >>> s = algebra.scalar(5)
        """
        # Import here to avoid circular imports at runtime
        from kingdon.multivector import MultiVector
        
        if isinstance(value, (list, tuple, np.ndarray)):
            if len(value) != 1:
                raise ValueError(f"Expected a single scalar value, got {len(value)} values")
            value = value[0]
        
        # Create a new MultiVector with scalar grade
        mv = object.__new__(MultiVector)
        mv.algebra = self
        mv._keys = (0,)
        mv._values = [value]
        mv.shape = ()
        
        return mv

    def vector(self, values: Union[List[Any], np.ndarray, Dict[int, Any]], 
               name: Optional[str] = None) -> 'MultiVector':
        """
        Create a vector multivector (grade 1).
        
        Args:
            values: Coefficients for the vector components
            name: Optional name for symbolic representation
            
        Returns:
            A new vector multivector
            
        Example:
            >>> # Create a vector with components [1, 2, 3]
            >>> v = algebra.vector([1, 2, 3])
        """
        # Import here to avoid circular imports at runtime
        from kingdon.multivector import MultiVector
        
        # Create a new MultiVector with vector grade
        mv = object.__new__(MultiVector)
        mv.algebra = self
        
        # Dimension determines the number of vector components
        d = self.d
        
        # Generate keys for grade-1 basis blades
        keys = tuple(1 << i for i in range(d))
        
        # Set values based on input type
        if isinstance(values, dict):
            # Dictionary input maps keys to values
            val_list = [0] * d
            for k, v in values.items():
                if 0 <= k < d:
                    val_list[k] = v
            mv._values = val_list
        else:
            # List or array input directly specifies values
            if not isinstance(values, (list, tuple, np.ndarray)):
                raise ValueError(f"Expected a list of vector components, got {type(values)}")
            if len(values) != d:
                raise ValueError(f"Expected {d} vector components, got {len(values)}")
            mv._values = list(values)
        
        mv._keys = keys
        mv.shape = ()
        
        return mv

    def bivector(self, values: Union[List[Any], np.ndarray, Dict[int, Any]], 
                 name: Optional[str] = None) -> 'MultiVector':
        """
        Create a bivector multivector (grade 2).
        
        Args:
            values: Coefficients for the bivector components
            name: Optional name for symbolic representation
            
        Returns:
            A new bivector multivector
            
        Example:
            >>> # Create a bivector
            >>> B = algebra.bivector([1, 0, 1])
        """
        # Import here to avoid circular imports at runtime
        from kingdon.multivector import MultiVector
        
        mv = object.__new__(MultiVector)
        mv.algebra = self
        
        d = self.d
        
        # Generate keys for grade-2 basis blades
        if d >= 3 and isinstance(values, (list, tuple, np.ndarray)) and len(values) == 3:
            # Special case for 3D: (e12, e13, e23) -> (3, 5, 6)
            keys = (3, 5, 6)
        else:
            # General case: all grade-2 basis blades
            keys = tuple(i for i in range(2**d) if bin(i).count('1') == 2)
        
        # Set values
        if isinstance(values, dict):
            val_list = [0] * len(keys)
            for i, k in enumerate(keys):
                val_list[i] = values.get(k, 0)
            mv._values = val_list
        else:
            if len(values) != len(keys):
                raise ValueError(f"Expected {len(keys)} bivector components, got {len(values)}")
            mv._values = list(values)
        
        mv._keys = keys
        mv.shape = ()
        
        return mv
    
    def trivector(self, values: Union[List[Any], np.ndarray, Dict[int, Any]], 
                 name: Optional[str] = None) -> 'MultiVector':
        """
        Create a trivector multivector (grade 3).
        
        Args:
            values: Coefficients for the trivector components
            name: Optional name for symbolic representation
            
        Returns:
            A new trivector multivector
            
        Example:
            >>> # Create a trivector
            >>> T = algebra.trivector([1])
        """
        return self.multivector(values=values, grades=(3,), name=name)

    def pseudoscalar(self, value: Any = 1, name: Optional[str] = None) -> 'MultiVector':
        """
        Create the pseudoscalar (highest grade element) for this algebra.
        
        Args:
            value: Coefficient for the pseudoscalar (default: 1)
            name: Optional name for symbolic representation
            
        Returns:
            The pseudoscalar multivector
            
        Example:
            >>> # Create the pseudoscalar
            >>> I = algebra.pseudoscalar()
        """
        # Clone the existing pseudoscalar and multiply by value
        if value == 1:
            return self.pss
        else:
            return self.pss * value

    def evenmv(self, values: Union[List[Any], np.ndarray], 
              name: Optional[str] = None) -> 'MultiVector':
        """
        Create a multivector with even grades only (0, 2, 4, ...).
        
        Args:
            values: Coefficients for the components
            name: Optional name for symbolic representation
            
        Returns:
            A new multivector with even grades
            
        Example:
            >>> # Create an even-grade multivector
            >>> R = algebra.evenmv([1, 0, 0, 1])
        """
        # Get all even grades up to the dimension
        even_grades = tuple(range(0, self.d + 1, 2))
        return self.multivector(values=values, grades=even_grades, name=name)

    def oddmv(self, values: Union[List[Any], np.ndarray], 
             name: Optional[str] = None) -> 'MultiVector':
        """
        Create a multivector with odd grades only (1, 3, 5, ...).
        
        Args:
            values: Coefficients for the components
            name: Optional name for symbolic representation
            
        Returns:
            A new multivector with odd grades
            
        Example:
            >>> # Create an odd-grade multivector
            >>> O = algebra.oddmv([1, 1, 0, 0])
        """
        # Get all odd grades up to the dimension
        odd_grades = tuple(range(1, self.d + 1, 2))
        return self.multivector(values=values, grades=odd_grades, name=name)

    def grade(self, mv: 'MultiVector', *grades: int) -> 'MultiVector':
        """
        Extract specified grades from a multivector.
        
        Args:
            mv: Input multivector
            *grades: Grades to extract (0=scalar, 1=vector, 2=bivector, etc.)
            
        Returns:
            A new multivector containing only the specified grades
            
        Example:
            >>> # Extract vector part from a multivector
            >>> v_part = algebra.grade(mv, 1)
            >>> # Extract scalar and bivector parts
            >>> sb_part = algebra.grade(mv, 0, 2)
        """
        return mv.grade(*grades)
    
    def blade(self, indices: List[int], value: Any = 1, name: Optional[str] = None) -> 'MultiVector':
        """
        Create a specific basis blade by indices.
        
        Args:
            indices: The indices of the basis vectors to wedge together
                (e.g., [1, 2] creates the e₁₂ blade)
            value: Coefficient for the blade (default: 1)
            name: Optional name for symbolic representation
            
        Returns:
            A multivector representing the specified basis blade
            
        Example:
            >>> # Create the e₁₃ blade with coefficient 2
            >>> B = algebra.blade([1, 3], 2)
        """
        if not indices:
            return self.scalar(value, name)
        
        # Try using the blades attribute
        if hasattr(self, 'blades'):
            # Create blade name (e.g., "e12" for e₁₂)
            sorted_indices = sorted(indices)
            blade_name = f"e{''.join(str(i) for i in sorted_indices)}"
            
            try:
                # Try attribute access first
                blade = getattr(self.blades, blade_name)
                return blade * value
            except AttributeError:
                try:
                    # Try dictionary access
                    blade = self.blades[blade_name]
                    return blade * value
                except (KeyError, TypeError):
                    pass
        
        # As a fallback, create blade using vector operations
        if len(indices) == 1:
            # For grade 1, just create a vector with the appropriate component
            idx = indices[0]
            vec_values = [0] * self.d
            adjusted_idx = idx - self.start_index
            if 0 <= adjusted_idx < self.d:
                vec_values[adjusted_idx] = value
            return self.vector(values=vec_values, name=name)
        else:
            # For higher grades, compose using outer product
            e_vectors = []
            for idx in indices:
                vec_values = [0] * self.d
                adjusted_idx = idx - self.start_index
                if 0 <= adjusted_idx < self.d:
                    vec_values[adjusted_idx] = 1
                e_vectors.append(self.vector(values=vec_values))
            
            # Use the outer product to combine the vectors
            result = e_vectors[0]
            for vec in e_vectors[1:]:
                result = result ^ vec
            
            return result * value

    def rotor(self, angle: float, plane_indices: Tuple[int, int], 
             name: Optional[str] = None) -> 'MultiVector':
        """
        Create a rotor representing a rotation in a specific plane.
        
        A rotor is of the form R = cos(θ/2) + sin(θ/2)Bᵢⱼ where Bᵢⱼ is a unit bivector.
        
        Args:
            angle: Rotation angle in radians
            plane_indices: Tuple of two indices defining the rotation plane
                (e.g., (1, 2) for rotation in the e₁e₂ plane)
            name: Optional name for symbolic representation
            
        Returns:
            A multivector representing the rotor
            
        Example:
            >>> import math
            >>> # Create a rotor for 90° rotation in the xy-plane
            >>> R = algebra.rotor(angle=math.pi/2, plane_indices=(1, 2))
        """
        import math
        # Create scalar part cos(θ/2)
        scalar = self.scalar(value=math.cos(angle/2))
        
        # Create bivector part sin(θ/2)Bᵢⱼ
        i, j = sorted(plane_indices)
        
        # Generate the bivector key 
        # e.g. for (1,2) we need e12 which corresponds to binary key 3
        binary_key = (1 << (i-1)) | (1 << (j-1))
        bivector = self.blade(indices=[i, j], value=math.sin(angle/2))
        
        # Import here to avoid circular imports at runtime
        from kingdon.multivector import MultiVector
        if isinstance(scalar, MultiVector) and isinstance(bivector, MultiVector):
            # Get the scalar and bivector values
            scalar_value = scalar._values[0]
            bivector_key = bivector._keys[0]
            bivector_value = bivector._values[0]
            
            # Create a new multivector with both components
            rotor = MultiVector.fromkeysvalues(
                self,
                keys=(0, bivector_key),
                values=[scalar_value, bivector_value]
            )
        else:
            # Fallback to addition if direct creation isn't possible
            rotor = scalar + bivector
        
        # Set name if needed
        if name and hasattr(rotor, 'name'):
            rotor.name = name
            
        return rotor

    def translator(self, direction: List[float], distance: float = 1.0, 
                  name: Optional[str] = None) -> 'MultiVector':
        """
        Create a translator in Projective Geometric Algebra (PGA).
        
        A translator is of the form T = 1 + (d/2)e₀ᵢvᵢ where e₀ᵢ is a null vector direction.
        This only works in algebras with at least one null dimension (r ≥ 1).
        
        Args:
            direction: Direction vector components (excluding null component)
            distance: Translation distance
            name: Optional name for symbolic representation
            
        Returns:
            A multivector representing the translator
            
        Example:
            >>> # Create a translator for 2 units in the x direction
            >>> T = algebra.translator(direction=[1, 0, 0], distance=2)
        """
        # Validate that the algebra has at least one null dimension
        if self.r < 1:
            raise ValueError("Translator requires an algebra with at least one null dimension (r ≥ 1)")
        
        # Create scalar part (1)
        translator = self.scalar(value=1)
        
        # Create bivector parts (d/2)e₀ᵢvᵢ
        for i, component in enumerate(direction, start=1):
            if component == 0:
                continue
                
            # Create bivector e₀ᵢ with coefficient (d/2)*vᵢ
            blade = self.blade(indices=[0, i], value=(distance/2) * component)
            translator = translator + blade
        
        # Set name if needed
        if name and hasattr(translator, 'name'):
            translator.name = name
            
        return translator

    def reflector(self, normal: List[float], name: Optional[str] = None) -> 'MultiVector':
        """
        Create a reflector (versor representing reflection in a hyperplane).
        
        A reflector is simply a normalized vector perpendicular to the reflection hyperplane.
        
        Args:
            normal: Normal vector components defining the reflection hyperplane
            name: Optional name for symbolic representation
            
        Returns:
            A multivector representing the reflector
            
        Example:
            >>> # Create a reflector for the xy-plane (normal = z-axis)
            >>> R = algebra.reflector(normal=[0, 0, 1])
        """
        import math
        
        # Validate input - but be more flexible with dimension
        if len(normal) < self.d - 1:
            raise ValueError(f"Normal vector must have at least {self.d-1} components, got {len(normal)}")
        
        # For PGA (r=1), we need to handle the special case where the normal might not include the homogeneous coordinate
        if self.r == 1 and len(normal) == self.d - 1:
            # Add a 0 for the homogeneous coordinate
            normal = list(normal) + [0]
        
        # Pad with zeros if needed
        if len(normal) < self.d:
            normal = list(normal) + [0] * (self.d - len(normal))
        
        # Check if the normal vector is zero
        norm_sq = sum(x*x for x in normal)
        if norm_sq == 0:
            raise ValueError("Normal vector cannot be zero")
        
        # Create the normal vector
        vector = self.vector(values=normal, name=name)
        
        # Calculate the norm manually for normalization
        norm = math.sqrt(norm_sq)
        
        # Import here to avoid circular imports at runtime
        from kingdon.multivector import MultiVector
        if isinstance(vector, MultiVector):
            # Create a new multivector with normalized values
            normalized_values = [v / norm for v in vector._values]
            return MultiVector.fromkeysvalues(self, vector._keys, normalized_values)
        else:
            # Fallback to division
            return vector / norm

    def motor(self, rotation: Optional[Tuple[float, Tuple[int, int]]] = None,
             translation: Optional[Tuple[List[float], float]] = None,
             name: Optional[str] = None) -> 'MultiVector':
        """
        Create a motor (combined rotation and translation) in PGA.
        
        A motor is the product of a rotor and a translator and represents a rigid body motion.
        
        Args:
            rotation: Optional tuple (angle, plane_indices) for rotation component
            translation: Optional tuple (direction, distance) for translation component
            name: Optional name for symbolic representation
            
        Returns:
            A multivector representing the motor
            
        Example:
            >>> import math
            >>> # Create a motor for 90° rotation in xy-plane followed by translation in x
            >>> M = algebra.motor(
            ...     rotation=(math.pi/2, (1, 2)),
            ...     translation=([1, 0, 0], 2)
            ... )
        """
        # Check if either rotation or translation is provided
        if not rotation and not translation:
            return self.scalar(value=1, name=name)
        
        # Check if the algebra has a null dimension for translation
        if translation and self.r < 1:
            raise ValueError("Translation requires an algebra with at least one null dimension (r ≥ 1)")
        
        # Start with identity motor
        motor = self.scalar(value=1)
        
        # Apply rotation if specified
        if rotation:
            angle, plane_indices = rotation
            rotor = self.rotor(angle, plane_indices)
            motor = motor * rotor
        
        # Apply translation if specified
        if translation:
            direction, distance = translation
            translator = self.translator(direction, distance)
            motor = motor * translator
        
        # Set name if needed
        if name and hasattr(motor, 'name'):
            motor.name = name
            
        return motor

    def dual(self, mv: 'MultiVector', kind: str = 'auto') -> 'MultiVector':
        """
        Compute the dual of a multivector.
        
        The dual maps a multivector to its orthogonal complement using the pseudoscalar.
        
        Args:
            mv: Input multivector
            kind: Type of dual to compute ('polarity', 'hodge', or 'auto')
            
        Returns:
            A new multivector representing the dual
            
        Example:
            >>> # Compute the dual of a vector
            >>> v_dual = algebra.dual(v)
        """
        return mv.dual(kind=kind)

    def undual(self, mv: 'MultiVector', kind: str = 'auto') -> 'MultiVector':
        """
        Compute the undual (inverse of dual) of a multivector.
        
        Args:
            mv: Input multivector
            kind: Type of undual to compute ('polarity', 'hodge', or 'auto')
            
        Returns:
            A new multivector representing the undual
            
        Example:
            >>> # Compute the undual of a bivector
            >>> b_undual = algebra.undual(b)
        """
        return mv.undual(kind=kind)

    @classmethod
    def from_signature(cls, signature: List[int]) -> 'Algebra':
        """
        Create an algebra from a custom signature.
        
        Args:
            signature: List of 1, -1, and 0 values defining the signature
            
        Returns:
            A new Algebra instance
            
        Example:
            >>> # Create a custom algebra with signature [1, 1, -1, 0]
            >>> alg = Algebra.from_signature([1, 1, -1, 0])
        """
        counts = Counter(signature)
        p, q, r = counts[1], counts[-1], counts[0]
        return cls(p=p, q=q, r=r, signature=signature)

    @classmethod
    def from_name(cls, name: str) -> 'Algebra':
        """
        Create a common algebra by name.
        
        Args:
            name: Name of the algebra to create
            
        Returns:
            A new Algebra instance
            
        Raises:
            ValueError: If the algebra name is not recognized
            
        Example:
            >>> # Create a 3D Conformal Geometric Algebra
            >>> cga3d = Algebra.from_name("CGA3D")
        """
        name = name.upper()
        
        algebras = {
            "CGA2D": (3, 1, 0),   # 2D Conformal Geometric Algebra
            "CGA3D": (4, 1, 0),   # 3D Conformal Geometric Algebra
            "PGA2D": (2, 0, 1),   # 2D Projective Geometric Algebra
            "PGA3D": (3, 0, 1),   # 3D Projective Geometric Algebra
            "VGA2D": (2, 0, 0),   # 2D Vector Geometric Algebra
            "VGA3D": (3, 0, 0),   # 3D Vector Geometric Algebra
            "STA": (1, 3, 0),     # Spacetime Algebra
            "APS": (0, 3, 0)      # Algebra of Physical Space
        }
        
        if name not in algebras:
            raise ValueError(f"Unknown algebra name: {name}")
        
        p, q, r = algebras[name]
        return cls(p=p, q=q, r=r)


# Convenience factory functions for common algebras
def PGA2D() -> Algebra:
    """Create a 2D Projective Geometric Algebra."""
    return Algebra(p=2, q=0, r=1)

def PGA3D() -> Algebra:
    """Create a 3D Projective Geometric Algebra."""
    return Algebra(p=3, q=0, r=1)

def CGA2D() -> Algebra:
    """Create a 2D Conformal Geometric Algebra."""
    return Algebra(p=3, q=1, r=0)

def CGA3D() -> Algebra:
    """Create a 3D Conformal Geometric Algebra."""
    return Algebra(p=4, q=1, r=0)

def VGA3D() -> Algebra:
    """Create a 3D Vector Geometric Algebra."""
    return Algebra(p=3, q=0, r=0)

def STA() -> Algebra:
    """Create a Spacetime Algebra."""
    return Algebra(p=1, q=3, r=0)

def APS() -> Algebra:
    """Create an Algebra of Physical Space."""
    return Algebra(p=0, q=3, r=0)
>>>>>>> 5d1776bf2c10cd778122088d9664260375c22a81
