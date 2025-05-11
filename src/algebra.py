# -*- coding: utf-8 -*-
"""
Unified Algebra Module for Geometric Algebra
(Docstring unchanged)
"""

from __future__ import annotations

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


from kingdon.operator_dict import OperatorDict, UnaryOperatorDict, BladeDict, Registry, AlgebraError
from kingdon.multivector   import MultiVector
from kingdon.codegen       import CodegenOutput

# Typing Imports
from typing import List, Dict, Optional, Any, Tuple, Callable, Union, Set, Type, TypeVar, cast, TYPE_CHECKING, ClassVar, Sequence # Added Sequence

# Third-party Imports
import numpy as np
import sympy

# Setup logging
log = logging.getLogger(__name__)



# Type variable for operator dictionaries
T = TypeVar('T', bound=OperatorDict)

# Field initializer for operator dictionaries
operation_field = partial(field, default_factory=dict, init=False, repr=False, compare=False)


@dataclass
class Algebra:
    """ Comprehensive Algebra class for Geometric (Clifford) Algebra. """
    p: int = 0
    q: int = 0
    r: int = 0
    signature: Optional[List[int]] = None
    start_index: Optional[int] = None
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
    registry: Dict[str, Union[OperatorDict, UnaryOperatorDict, Callable]] = field(init=False, repr=False, compare=False, default_factory=dict) # Allow Callables for register
    numspace: Dict[str, Any] = field(init=False, repr=False, compare=False, default_factory=dict)

    # Operator dictionaries (Definitions unchanged)
    gp: OperatorDict = operation_field()
    sw: OperatorDict = operation_field()
    cp: OperatorDict = operation_field()
    acp: OperatorDict = operation_field()
    ip: OperatorDict = operation_field()
    sp: OperatorDict = operation_field()
    lc: OperatorDict = operation_field()
    rc: OperatorDict = operation_field()
    op: OperatorDict = operation_field()
    rp: OperatorDict = operation_field()
    proj: OperatorDict = operation_field()
    add: OperatorDict = operation_field()
    sub: OperatorDict = operation_field()
    div: OperatorDict = operation_field()
    inv: UnaryOperatorDict = operation_field()
    neg: UnaryOperatorDict = operation_field()
    reverse: UnaryOperatorDict = operation_field()
    involute: UnaryOperatorDict = operation_field()
    conjugate: UnaryOperatorDict = operation_field()
    sqrt: UnaryOperatorDict = operation_field()
    polarity: UnaryOperatorDict = operation_field()
    unpolarity: UnaryOperatorDict = operation_field()
    hodge: UnaryOperatorDict = operation_field()
    unhodge: UnaryOperatorDict = operation_field()
    normsq: UnaryOperatorDict = operation_field()
    outerexp: UnaryOperatorDict = operation_field()
    outersin: UnaryOperatorDict = operation_field()
    outercos: UnaryOperatorDict = operation_field()
    outertan: UnaryOperatorDict = operation_field()

    # --- Metadata for Operators (Unchanged) ---
    _operator_metadata: ClassVar[Dict[str, Dict[str, Any]]] = {
        'gp': {'codegen': 'codegen_gp', 'cls': OperatorDict}, 'sw': {'codegen': 'codegen_sw', 'cls': OperatorDict},
        'cp': {'codegen': 'codegen_cp', 'cls': OperatorDict}, 'acp': {'codegen': 'codegen_acp', 'cls': OperatorDict},
        'ip': {'codegen': 'codegen_ip', 'cls': OperatorDict}, 'sp': {'codegen': 'codegen_sp', 'cls': OperatorDict},
        'lc': {'codegen': 'codegen_lc', 'cls': OperatorDict}, 'rc': {'codegen': 'codegen_rc', 'cls': OperatorDict},
        'op': {'codegen': 'codegen_op', 'cls': OperatorDict}, 'rp': {'codegen': 'codegen_rp', 'cls': OperatorDict},
        'proj': {'codegen': 'codegen_proj', 'cls': OperatorDict}, 'add': {'codegen': 'codegen_add', 'cls': OperatorDict},
        'sub': {'codegen': 'codegen_sub', 'cls': OperatorDict}, 'div': {'codegen': 'codegen_div', 'cls': OperatorDict},
        'inv': {'codegen': 'codegen_inv', 'cls': UnaryOperatorDict}, 'neg': {'codegen': 'codegen_neg', 'cls': UnaryOperatorDict},
        'reverse': {'codegen': 'codegen_reverse', 'cls': UnaryOperatorDict}, 'involute': {'codegen': 'codegen_involute', 'cls': UnaryOperatorDict},
        'conjugate': {'codegen': 'codegen_conjugate', 'cls': UnaryOperatorDict}, 'sqrt': {'codegen': 'codegen_sqrt', 'cls': UnaryOperatorDict},
        'polarity': {'codegen': 'codegen_polarity', 'cls': UnaryOperatorDict}, 'unpolarity': {'codegen': 'codegen_unpolarity', 'cls': UnaryOperatorDict},
        'hodge': {'codegen': 'codegen_hodge', 'cls': UnaryOperatorDict}, 'unhodge': {'codegen': 'codegen_unhodge', 'cls': UnaryOperatorDict},
        'normsq': {'codegen': 'codegen_normsq', 'cls': UnaryOperatorDict}, 'outerexp': {'codegen': 'codegen_outerexp', 'cls': UnaryOperatorDict},
        'outersin': {'codegen': 'codegen_outersin', 'cls': UnaryOperatorDict}, 'outercos': {'codegen': 'codegen_outercos', 'cls': UnaryOperatorDict},
        'outertan': {'codegen': 'codegen_outertan', 'cls': UnaryOperatorDict},
    }

    # Add _grade_index_map caching to __post_init__ for efficiency
    def __post_init__(self) -> None:
        # Imports moved here to avoid potential circularity at module level
        # Need to handle potential import errors if modules aren't found
        try:
            from kingdon.multivector import MultiVector
            from kingdon.operator_dict import BladeDict
        except ImportError as e:
             log.error(f"Failed to import MultiVector or BladeDict in Algebra.__post_init__: {e}")
             # Depending on severity, raise error or try to continue with limited functionality
             raise ImportError(f"Kingdon core components (MultiVector/BladeDict) not found: {e}") from e

        self.numspace = {'np': np, 'MultiVector': MultiVector, 'cos': np.cos, 'sin': np.sin, 'sqrt': np.sqrt,
                         'math': math} # Added math
        self._setup_signature()
        self._setup_basis_mappings()
        self._precompute_grade_indices()
        self._compute_multiplication_signs()

        # --- Add caching for sort key map ---
        self._key_to_grade_index_map = {}
        # Use range(len(self)) to iterate through 0 to 2^d - 1
        # Groupby requires sorted input
        all_keys_iterable = range(len(self))
        key_func = lambda k: bin(k).count('1')
        all_keys_sorted_by_grade = sorted(all_keys_iterable, key=key_func)

        for grade, group in groupby(all_keys_sorted_by_grade, key=key_func):
             grade_keys = list(group)
             for index, key in enumerate(grade_keys):
                 # Store grade and index within grade (as float for safety)
                 self._key_to_grade_index_map[key] = (grade, float(index))
        # --- End caching ---

        self._initialize_blade_dictionary(BladeDict, MultiVector) # Sets self.pss
        self._create_registry()
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
        if self.start_index == 0: self.canon2bin["e0"] = 0 # Ensure e0 mapping if starting at 0

    def _blade2canon(self, basis_blade: str) -> Tuple[str, int]:
        if basis_blade == "e" or (self.start_index == 0 and basis_blade == "e0"): return "e", 0
        match_default = re.match(r'^e([0-9]+)$', basis_blade)
        if match_default:
            indices_str = match_default.group(1)
            try: indices = [int(c) for c in indices_str]
            except ValueError as e: raise ValueError(f"Invalid index character in blade string: {basis_blade}. {e}") from e
            if len(set(indices)) != len(indices): raise ValueError(f"Duplicate indices in blade name: {basis_blade}")
            min_idx, max_idx = self.start_index, self.start_index + self.d - 1
            if not all(min_idx <= i <= max_idx for i in indices): raise ValueError(f"Indices in '{basis_blade}' must be between {min_idx} and {max_idx}")
            sorted_indices = sorted(indices); canonical_name = f"e{''.join(map(str, sorted_indices))}"
            # Check if canonical name exists before calculating swaps (optimization/validation)
            if canonical_name not in self.canon2bin and canonical_name != 'e':
                 # This case implies indices formed a valid combination, but it wasn't generated? Should not happen.
                 raise ValueError(f"Internal error: Canonical name '{canonical_name}' for indices {sorted_indices} not found in algebra.")
            # Calculate swaps using bubble sort logic (inefficient but simple)
            current_indices = list(indices); swaps = 0; n = len(current_indices)
            for i in range(n):
                for j in range(n - 1 - i):
                    if current_indices[j] > current_indices[j+1]: current_indices[j], current_indices[j+1] = current_indices[j+1], current_indices[j]; swaps += 1
            return canonical_name, swaps
        elif self.basis_names:
            matched_names = []
            remaining_blade = basis_blade
            # Sort custom basis names by length descending to match longest first
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
            # Sort matched custom names based on their original index order
            name_to_index = {name: i for i, name in enumerate(self.basis_names)}
            sorted_matched_names = sorted(matched_names, key=lambda name: name_to_index[name])
            # Construct canonical name by joining sorted custom names
            # NOTE: This assumes concatenation is the desired canonical form for custom names.
            # A more robust system might map combined custom names to a standard 'e...' form.
            # For now, we map back to the 'e<indices>' form based on the keys.
            try:
                 combined_key = reduce(operator.xor, (self.canon2bin[n] for n in sorted_matched_names))
                 canonical_name = self.bin2canon.get(combined_key)
                 if canonical_name is None:
                     # Fallback if direct key lookup fails (should not happen if setup is correct)
                     log.warning(f"Could not find canonical name for key {combined_key} from custom names {sorted_matched_names}. Using key directly.")
                     canonical_name = f"key_{combined_key}" # Placeholder
            except KeyError as e:
                raise ValueError(f"Internal error: Basis name '{e}' used in '{basis_blade}' not found in algebra's canon2bin mapping.")

            # Calculate swaps
            current_names = list(matched_names)
            swaps = 0; n = len(current_names)
            for i in range(n):
                for j in range(n - 1 - i):
                    if name_to_index[current_names[j]] > name_to_index[current_names[j+1]]:
                        current_names[j], current_names[j+1] = current_names[j+1], current_names[j]
                        swaps += 1
            return canonical_name, swaps
        else:
             # Check if it's a pre-existing canonical name (e.g., "e12")
             if basis_blade in self.canon2bin: return basis_blade, 0
             # Otherwise, invalid format
             raise ValueError(f"Invalid basis blade format: '{basis_blade}'. Use 'e' + indices or provide custom basis_names.")

    def _precompute_grade_indices(self) -> None:
        self.indices_for_grades = defaultdict(list)
        # Use range(len(self)) which is 0 to 2^d - 1
        all_keys_iterable = range(len(self))
        key_func = lambda k: bin(k).count('1')
        # Groupby requires sorted input
        all_keys_sorted = sorted(all_keys_iterable, key=key_func)
        for grade, group in groupby(all_keys_sorted, key=key_func):
            # Store as tuple key (grade,)
            self.indices_for_grades[(grade,)] = list(group)

    def _compute_multiplication_signs(self) -> None:
        self.signs = {};
        algebra_len = len(self) # 2**d
        for i in range(algebra_len):
            for j in range(algebra_len):
                if i == 0 or j == 0: self.signs[(i, j)] = 1; continue
                i_indices = {k for k in range(self.d) if (i >> k) & 1}
                j_indices = {k for k in range(self.d) if (j >> k) & 1}
                common_indices = i_indices.intersection(j_indices)
                i_only_indices = i_indices.difference(common_indices)
                j_only_indices = j_indices.difference(common_indices)
                # Calculate swaps needed to bring common elements together and sort remaining
                swaps = 0
                for i_idx in i_only_indices:
                    swaps += sum(1 for j_idx in j_only_indices if i_idx > j_idx) # Count j's that should be before i

                sign = (-1)**swaps
                for k in common_indices:
                     if k >= len(self.signature): # Check bounds
                          log.error(f"Index {k} out of bounds for signature length {len(self.signature)} when calculating sign for ({i}, {j}).")
                          sign = 0; break # Treat as error / zero product
                     sig_val = self.signature[k]
                     if sig_val == 0: sign = 0; break
                     sign *= sig_val
                # Store only non-zero results, can infer zero otherwise
                if sign != 0: self.signs[(i, j)] = int(sign)

    def _initialize_blade_dictionary(self, BladeDictClass: Type['BladeDict'], MVClass: Type['MultiVector']) -> None:
        # Lazy initialization heuristic
        lazy_init = hasattr(self, 'd') and self.d > 6
        self.blades = BladeDictClass(self, lazy=lazy_init)
        # Initialize Pseudoscalar (pss)
        if hasattr(self, 'd') and self.d >= 0:
            if self.d == 0: # Special case for d=0 algebra (scalars only)
                 pseudoscalar_key = 0
                 pseudoscalar_name = "e"
            else:
                 pseudoscalar_key = (1 << self.d) - 1
                 pseudoscalar_name = self.bin2canon.get(pseudoscalar_key)

            if pseudoscalar_name:
                 # Access via BladeDict __getitem__ to trigger creation if lazy
                 try:
                      self.pss = self.blades[pseudoscalar_name]
                 except Exception as e:
                      log.error(f"Failed to get/create pss '{pseudoscalar_name}' via BladeDict: {e}. Creating manually.")
                      # Manually create if BladeDict access fails
                      self.pss = MVClass.fromkeysvalues(self, keys=(pseudoscalar_key,), values=[1])
            else:
                log.warning(f"Could not find canonical name for pseudoscalar key {pseudoscalar_key}. Creating manually.")
                self.pss = MVClass.fromkeysvalues(self, keys=(pseudoscalar_key,), values=[1])
        else:
             # Handle case where dimension d is not set (should not happen after __post_init__)
             log.error("Algebra dimension 'd' not set before initializing pss.")
             # Set pss to None or raise error? Let's set to None.
             self.pss = None # Or raise AttributeError("Algebra not properly initialized.")


    def _create_registry(self) -> None:
        # Need to import these locally if not already available
        try:
             from kingdon.operator_dict import OperatorDict, UnaryOperatorDict
             import kingdon.codegen as codegen_module
        except ImportError as e:
             log.error(f"Failed to import OperatorDict/UnaryOperatorDict/codegen in _create_registry: {e}")
             self.registry = {} # Initialize empty registry on failure
             return # Cannot proceed

        self.registry = {}
        for op_name, meta in self._operator_metadata.items():
            codegen_func_name = meta.get('codegen'); operator_cls_name = meta.get('cls', OperatorDict)
            if codegen_func_name:
                codegen_func = getattr(codegen_module, codegen_func_name, None)
                if codegen_func is None:
                     log.warning(f"Codegen function '{codegen_func_name}' not found for operator '{op_name}'. Skipping registration.")
                     continue
                try:
                     operator_cls = UnaryOperatorDict if operator_cls_name == UnaryOperatorDict else OperatorDict
                     operator_instance = operator_cls(name=op_name, algebra=self, codegen=codegen_func)
                     self.registry[op_name] = operator_instance
                     setattr(self, op_name, operator_instance)
                except Exception as e_create:
                     log.error(f"Failed to create/register operator '{op_name}': {e_create}")


    def _compile_gp_func(self) -> Callable[[int, int], Tuple[int, int]]:
        @lru_cache(maxsize=None)
        def gp_basis(b1_idx: int, b2_idx: int) -> Tuple[int, int]:
            """ Computes basis_out_idx, sign = gp(basis[b1_idx], basis[b2_idx])."""
            sign = self.signs.get((b1_idx, b2_idx), 0)
            # Resulting blade index is XOR of input indices
            result_idx = b1_idx ^ b2_idx
            return result_idx, sign
        return gp_basis

    # ========================================================================
    # Helper Methods
    # ========================================================================

    def _parse_composite_custom_key(self, key_str: str) -> Optional[List[int]]:
        # (Unchanged from previous state - relies on self.basis_names, self.canon2bin)
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
        # Return list of integer keys corresponding to matched custom names
        return matched_keys

    # --- ADDED: Definition for _canonical_key_sort_func ---
    def _canonical_key_sort_func(self) -> Callable[[int], Tuple[int, float]]:
        """
        Returns a sort key function for canonical ordering (grade, then index within grade).
        Uses a precomputed map for efficiency.
        """
        # Ensure the map is built (should be done in __post_init__)
        if not hasattr(self, '_key_to_grade_index_map'):
             # Fallback build if somehow missed in post_init (should not happen)
             print("Warning: _key_to_grade_index_map not found, rebuilding...")
             # (Rebuild logic - copy from __post_init__)
             self._key_to_grade_index_map = {}
             all_keys_iterable = range(len(self))
             key_func = lambda k: bin(k).count('1')
             all_keys_sorted_by_grade = sorted(all_keys_iterable, key=key_func)
             for grade, group in groupby(all_keys_sorted_by_grade, key=key_func):
                  grade_keys = list(group)
                  for index, key in enumerate(grade_keys):
                      self._key_to_grade_index_map[key] = (grade, float(index))

        # Use the map built during __post_init__
        _map = self._key_to_grade_index_map # Local reference for the returned function

        def sort_key(k: int) -> Tuple[int, float]:
            """Sorts by grade, then index within that grade using precomputed map."""
            # Provide a fallback: sort by grade, then numerical key value if not found
            return _map.get(k, (bin(k).count('1'), float(k)))

        return sort_key
    # --- END ADDED ---

    # ========================================================================
    # Multivector Creation Methods (Using _process_creation_args)
    # ========================================================================

    def _process_creation_args(self, *args: Any, **kwargs: Any) -> Tuple[Optional[Tuple[int, ...]], Optional[Union[List[Any], np.ndarray, Dict[str, Any]]], Optional[Tuple[int, ...]], Optional[str], Optional[Any]]:
        """Internal helper to parse arguments for multivector creation."""
        # Pop arguments from kwargs first
        values = kwargs.pop('values', None)
        keys = kwargs.pop('keys', None)
        grades = kwargs.pop('grades', None)
        name = kwargs.pop('name', None)
        symbolcls = kwargs.pop('symbolcls', None)

        # --- Argument Conflict Checks ---
        # (Combined checks for clarity)
        has_keys = keys is not None
        has_grades = grades is not None
        has_kwargs_blades = bool(kwargs) # Check if remaining kwargs exist (implies blade keys)
        has_name = name is not None
        has_values = values is not None

        if has_keys and has_grades: raise ValueError("Cannot specify both 'keys' and 'grades' arguments.")
        if has_keys and has_kwargs_blades: raise ValueError("Cannot specify both 'keys' and keyword arguments for blades.")
        if has_grades and has_kwargs_blades: raise ValueError("Cannot specify 'grades' and keyword arguments for blades.")
        if has_name and has_kwargs_blades: raise ValueError("Cannot specify both 'name' (symbolic) and keyword arguments for blades.")
        if has_name and has_values: raise ValueError("Cannot specify both 'name' (symbolic) and explicit 'values'.")

        # --- Process Positional Args ---
        if args:
            if len(args) > 1: raise TypeError(f"Too many positional arguments ({len(args)}), expected 0 or 1.")
            arg = args[0]
            # Assign positional arg based on what's not already set by keywords
            if isinstance(arg, str) and not has_name: name = arg
            elif not has_values: values = arg
            elif isinstance(arg, (int, float, complex, sympy.Expr)) and not has_values and not has_name: values = arg # Allow scalar value
            else: raise TypeError(f"Unexpected positional argument type '{type(arg).__name__}' or argument already specified via keyword.")

        # --- Process Keyword Args (Blades) ---
        if has_kwargs_blades:
             # Ensure basis mappings are ready
             if not hasattr(self, 'canon2bin') or not self.canon2bin: self._setup_basis_mappings()
             blade_kwargs_processed = {}
             for k_str, v in kwargs.items():
                  # Validate blade name (standard 'e...', custom, or composite custom)
                  is_valid_blade = False
                  if isinstance(k_str, str):
                       try:
                            # Check if it's a valid canonical or custom name via _blade2canon
                            self._blade2canon(k_str) # This validates format/indices
                            is_valid_blade = True
                       except ValueError:
                            # If not standard, try parsing as composite custom key
                            if self._parse_composite_custom_key(k_str): is_valid_blade = True
                  elif k_str == 'e': is_valid_blade = True # Explicitly allow 'e'

                  if is_valid_blade:
                      blade_kwargs_processed[k_str] = v
                  else:
                      # Raise error for unexpected keyword args not matching blades
                      raise TypeError(f"{type(self).__name__}.multivector() got an unexpected keyword argument '{k_str}'. If specifying blades, use valid names.")

             if blade_kwargs_processed:
                 # Merge blade kwargs into the 'values' dict
                 if has_values:
                      if not isinstance(values, dict):
                           raise TypeError("Cannot specify both 'values' (as list/array/scalar) and keyword arguments for blades. Use a dictionary for 'values' or omit it.")
                      # Check for overlaps
                      overlapping_keys = set(values.keys()) & set(blade_kwargs_processed.keys())
                      if overlapping_keys: raise ValueError(f"Blade keys specified in both 'values' dict and keyword arguments: {overlapping_keys}")
                      values.update(blade_kwargs_processed)
                 else:
                      values = blade_kwargs_processed # Assign if values wasn't provided

        # Return processed/consolidated arguments
        return keys, values, grades, name, symbolcls


    def multivector(self, *args: Any, **kwargs: Any) -> 'MultiVector':
        """ Primary factory method to create a multivector in this algebra. """
        # Import here to avoid circularity at module level if needed
        try:
            from kingdon.multivector import MultiVector
        except ImportError as e:
             log.error(f"Failed to import MultiVector in Algebra.multivector: {e}")
             raise

        keys, values, grades, name, symbolcls_arg = self._process_creation_args(*args, **kwargs)
        symbolcls_arg = symbolcls_arg or self.codegen_symbolcls # Use algebra default if not specified

        # --- Variable Initialization ---
        final_keys: Optional[List[int]] = None # Use list initially for easier manipulation
        final_values: Union[List[Any], np.ndarray]
        target_len: int = -1 # Expected length based on keys/grades

        # --- Determine Final Keys and Target Length ---
        max_key_val = (1 << self.d) - 1 if self.d >= 0 else -1 # Maximum valid integer key

        if keys is not None: # Case 1: Explicit keys provided
            processed_keys_set = set() # Use set to handle duplicates gracefully
            keys_input = (keys,) if isinstance(keys, (int, str)) else keys # Normalize input

            if not isinstance(keys_input, (tuple, list, set)): # Allow set input as well
                raise TypeError(f"`keys` argument must be an int, string, or an iterable (tuple/list/set). Got {type(keys_input).__name__}.")

            for k in keys_input:
                if isinstance(k, int):
                    if not (0 <= k <= max_key_val): raise ValueError(f"Integer key {k} out of range [0, {max_key_val}] for dimension {self.d}.")
                    processed_keys_set.add(k)
                elif isinstance(k, str):
                    try:
                        canon_name, _ = self._blade2canon(k) # Validates format/indices
                        int_key = self.canon2bin.get(canon_name)
                        if int_key is None: raise KeyError(f"Internal error mapping canonical name '{canon_name}'.")
                        processed_keys_set.add(int_key)
                    except (ValueError, KeyError) as e: raise ValueError(f"Invalid basis blade name '{k}' provided in 'keys': {e}") from e
                else: raise TypeError(f"Items in `keys` must be integers or strings. Got type {type(k).__name__}.")

            # Convert set to sorted list using canonical sort func
            final_keys = sorted(list(processed_keys_set), key=self._canonical_key_sort_func())
            target_len = len(final_keys)
            # Handle d=0 case where keys might be just {0}
            if not final_keys and 0 in processed_keys_set and self.d == 0: final_keys = [0]; target_len = 1

        elif grades is not None: # Case 2: Grades provided
            grades_input = (grades,) if isinstance(grades, int) else grades # Normalize
            if not isinstance(grades_input, (tuple, list, set)): # Allow set
                raise TypeError(f"`grades` must be an int or an iterable (tuple/list/set). Got {type(grades_input).__name__}.")

            valid_grades = set()
            for g in grades_input:
                 if not isinstance(g, int) or not (0 <= g <= self.d): raise ValueError(f"Grade {g} invalid for dim {self.d}.")
                 valid_grades.add(g)

            collected_keys_set = set()
            for g in sorted(list(valid_grades)):
                grade_keys = self.indices_for_grades.get((g,), [])
                collected_keys_set.update(grade_keys)

            final_keys = sorted(list(collected_keys_set), key=self._canonical_key_sort_func())
            target_len = len(final_keys)
            if not final_keys and 0 in valid_grades and self.d == 0: final_keys = [0]; target_len = 1

        elif isinstance(values, dict): # Case 3: Values as Dict (already merged with kwargs)
             processed_keys_dict = {} # Map int_key -> original_key_input
             keys_to_process_set = set()
             for k_in, v in values.items():
                 int_key: Optional[int] = None
                 original_key_repr = repr(k_in) # Store representation for error messages
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
                     except (ValueError, KeyError) as e: raise KeyError(f"Invalid basis blade key '{k_in}': {e}") from e
                 elif isinstance(k_in, int):
                     if not (0 <= k_in <= max_key_val): raise KeyError(f"Integer key {k_in} out of range [0, {max_key_val}].")
                     int_key = k_in
                 else: raise TypeError(f"Invalid key type {type(k_in).__name__} in values dict. Use int or str.")

                 if int_key is None: raise RuntimeError(f"Internal error: Failed to determine integer key for input key '{k_in}'.")
                 if int_key in keys_to_process_set: raise ValueError(f"Duplicate key specified: '{k_in}' resolves to same component as input '{processed_keys_dict[int_key]}'.")
                 keys_to_process_set.add(int_key)
                 processed_keys_dict[int_key] = k_in # Store original key for value lookup later

             final_keys = sorted(list(keys_to_process_set), key=self._canonical_key_sort_func())
             target_len = len(final_keys)
             # Handle empty dict case - results in zero scalar
             if not final_keys: final_keys = [0]; target_len = 1

        elif name is not None: # Case 4: Symbolic multivector requested
             # Check if grades were also specified implicitly (should be caught by conflict checks)
             if grades is not None: pass # final_keys/target_len already set by Case 2
             else: # Full symbolic MV (all canonical blades)
                 # Use range(len(self)) to get all keys 0 to 2^d-1
                 final_keys = sorted(list(range(len(self))), key=self._canonical_key_sort_func())
                 target_len = len(self) # Length of full algebra

        else: # Case 5: Values as list/array/scalar or None (default)
            # Default to scalar 0 if values is None
            if values is None: values = 0

            if isinstance(values, (int, float, complex, sympy.Expr)):
                 # Interpret as scalar value if keys/grades not specified
                 final_keys = [0]; target_len = 1
            elif isinstance(values, (list, tuple, np.ndarray)): # Sequence input
                 val_len = len(values)
                 algebra_len = len(self)
                 if val_len == algebra_len: # Assume full canonical sequence
                     final_keys = sorted(list(range(algebra_len)), key=self._canonical_key_sort_func())
                     target_len = algebra_len
                 else: # Ambiguous length - require explicit keys/grades
                      raise ValueError(f"Length of sequence 'values' ({val_len}) is ambiguous. Provide explicit 'keys' or 'grades', or ensure length matches algebra size ({algebra_len}).")
            else: # Unhandled type for values when keys/grades not specified
                 raise TypeError(f"Unsupported type for 'values': {type(values).__name__}. Provide explicit 'keys' or 'grades', or use dict/sequence/scalar.")

        # --- Final Key/Length Consistency Check ---
        if final_keys is None or target_len == -1: raise RuntimeError("Internal state error: final_keys or target_len not determined.")
        # Ensure final_keys is a list before proceeding
        if not isinstance(final_keys, list): final_keys = list(final_keys)

        # --- Determine Final Values ---
        is_sympy_algebra = self.codegen_symbolcls != float # Check if algebra is symbolic
        default_zero = sympy.Integer(0) if is_sympy_algebra else 0.0

        if name is not None: # Case 4: Symbolic
             sym_func = symbolcls_arg or sympy.Symbol
             final_values_list = []
             for k in final_keys:
                 # Use bin2canon for suffix, provide fallback
                 blade_suffix = self.bin2canon.get(k, f'k{k}')
                 # Handle scalar case name
                 sym_name = name if k == 0 else f"{name}_{blade_suffix}"
                 try:
                     # Create symbol, assume scalar for key 0 if name only provided?
                     # Or require specific dict value for scalar? Assume sym for k=0 too.
                     final_values_list.append(sym_func(sym_name))
                 except TypeError as e_sym: raise TypeError(f"Failed to create symbol '{sym_name}' using symbol class {sym_func}: {e_sym}")
             final_values = final_values_list

        elif isinstance(values, dict): # Case 3: Dict/Keywords
             # Need mapping from int_key back to original input key for lookup
             int_to_orig_key = {k: v for v, k in processed_keys_dict.items()}
             final_values_map = {}
             for k_int in final_keys:
                  orig_key = int_to_orig_key.get(k_int)
                  if orig_key is None: # Should not happen if final_keys derived correctly
                       log.error(f"Internal error: Cannot find original key for integer key {k_int} in dict case.")
                       val_to_use = default_zero
                  else:
                       raw_val = values.get(orig_key, default_zero)
                       # Sympify only if symbolic algebra context, else keep as is
                       val_to_use = sympify(raw_val) if is_sympy_algebra else raw_val
                  final_values_map[k_int] = val_to_use
             # Handle the zero scalar case where final_keys might be [0] but values dict was empty
             if not values and final_keys == [0]: final_values = [default_zero]
             else: final_values = [final_values_map.get(k, default_zero) for k in final_keys]

        elif isinstance(values, (list, tuple, np.ndarray)): # Case 1, 2, 5 (Sequence)
             if len(values) != target_len: # Should be caught earlier, but double check
                 raise ValueError(f"Internal value length mismatch: expected {target_len}, got {len(values)}.")
             # Convert to list/ndarray, sympify elements if symbolic algebra
             processed_vals = [sympify(v) if is_sympy_algebra else v for v in values]
             final_values = np.array(processed_vals) if isinstance(values, np.ndarray) else list(processed_vals)

        elif isinstance(values, (int, float, complex, sympy.Expr)): # Case 5 (Scalar)
             if final_keys == [0] and target_len == 1:
                 # Sympify if symbolic algebra
                 final_values = [sympify(values) if is_sympy_algebra else values]
             else: # Should be caught earlier
                 raise ValueError(f"Scalar value '{values}' provided, but target keys {final_keys} indicate non-scalar structure.")
        else: # Should not be reached if logic above is correct
             raise TypeError(f"Internal error: Unhandled type for 'values': {type(values).__name__}")

        # --- Final Construction using fromkeysvalues ---
        # Ensure keys is tuple, values is list or ndarray
        final_keys_tuple = tuple(final_keys)
        return MultiVector.fromkeysvalues(self, keys=final_keys_tuple, values=final_values)

    # ========================================================================
    # Convenience Factory Methods (Calling self.multivector)
    # ========================================================================
    # (Methods scalar, vector, bivector, trivector, purevector, evenmv, oddmv
    #  blade, pseudoscalar, rotor, translator, reflector, motor remain unchanged
    #  from previous step, as they already call self.multivector)
    def scalar(self, value: Any = 0, *, name: Optional[str] = None) -> 'MultiVector':
        """Creates a scalar (grade 0) multivector."""
        # Pass value directly, multivector handles scalar case
        return self.multivector(values=value, grades=(0,), name=name)

    def vector(self, values: Union[Sequence[Any], np.ndarray, Dict[Union[str, int], Any]], *, name: Optional[str] = None) -> 'MultiVector':
        """Creates a vector (grade 1) multivector."""
        return self.multivector(values=values, grades=(1,), name=name)

    def bivector(self, values: Union[Sequence[Any], np.ndarray, Dict[Union[str, int], Any]], *, name: Optional[str] = None) -> 'MultiVector':
        """Creates a bivector (grade 2) multivector."""
        return self.multivector(values=values, grades=(2,), name=name)

    def trivector(self, values: Union[Sequence[Any], np.ndarray, Dict[Union[str, int], Any]], *, name: Optional[str] = None) -> 'MultiVector':
        """Creates a trivector (grade 3) multivector."""
        return self.multivector(values=values, grades=(3,), name=name)

    def purevector(self, values: Union[Sequence[Any], np.ndarray, Dict[Union[str, int], Any]], grade: int, *, name: Optional[str] = None) -> 'MultiVector':
        """Creates a pure grade multivector."""
        if not isinstance(grade, int) or not (0 <= grade <= self.d):
             raise ValueError(f"Grade must be an integer between 0 and {self.d}")
        return self.multivector(values=values, grades=(grade,), name=name)

    def evenmv(self, values: Union[Sequence[Any], np.ndarray, Dict[Union[str, int], Any]], *, name: Optional[str] = None) -> 'MultiVector':
        """Creates a multivector containing only even grades."""
        even_grades = tuple(g for g in range(0, self.d + 1) if g % 2 == 0)
        if not even_grades and self.d == 0: even_grades = (0,) # Handle d=0 case
        # Return zero scalar if no even grades exist (e.g., d=-1?)
        # elif not even_grades: return self.scalar(0, name=name)
        return self.multivector(values=values, grades=even_grades, name=name)

    def oddmv(self, values: Union[Sequence[Any], np.ndarray, Dict[Union[str, int], Any]], *, name: Optional[str] = None) -> 'MultiVector':
        """Creates a multivector containing only odd grades."""
        odd_grades = tuple(g for g in range(1, self.d + 1) if g % 2 != 0)
        # Return zero scalar if no odd grades exist (e.g., d=0)
        # if not odd_grades: return self.scalar(0, name=name)
        return self.multivector(values=values, grades=odd_grades, name=name)

    def blade(self, spec: Union[str, Sequence[int]], value: Any = 1, *, name: Optional[str] = None) -> 'MultiVector':
        """Creates a multivector representing a single basis blade."""
        # Use dict format for values to pass to self.multivector
        blade_values = {}
        sign = 1

        if isinstance(spec, (list, tuple, set)):
            indices = list(spec)
            if not indices: key_str = "e" # Scalar
            else:
                # Validate indices
                min_idx, max_idx = self.start_index, self.start_index + self.d - 1
                if not all(isinstance(i, int) and min_idx <= i <= max_idx for i in indices):
                     raise ValueError(f"Indices {indices} must be integers between {min_idx} and {max_idx}")
                if len(set(indices)) != len(indices): raise ValueError(f"Indices in {indices} must be unique.")
                # Get canonical name and sign factor
                sorted_indices = sorted(indices)
                key_str = f"e{''.join(map(str, sorted_indices))}"
                current_indices = list(indices); swaps = 0; n = len(current_indices)
                for i in range(n): # Bubble sort swap counting
                     for j in range(n - 1 - i):
                          if current_indices[j] > current_indices[j+1]: current_indices[j], current_indices[j+1] = current_indices[j+1], current_indices[j]; swaps += 1
                sign = (-1)**swaps
        elif isinstance(spec, str):
            try:
                canonical_name, swaps = self._blade2canon(spec)
                key_str = canonical_name
                sign = (-1)**swaps
            except ValueError as e: raise ValueError(f"Invalid blade specification string '{spec}': {e}")
        else: raise TypeError("Blade specification must be a string or an iterable of indices.")

        final_value = sign * value
        blade_values[key_str] = final_value
        # Pass dict to multivector factory
        try: return self.multivector(values=blade_values, name=name)
        except KeyError as e: raise ValueError(f"Could not create blade for spec '{spec}' (resolved to '{key_str}'): {e}")

    def pseudoscalar(self, value: Any = 1, *, name: Optional[str] = None) -> 'MultiVector':
        """Creates the pseudoscalar multivector, optionally scaled."""
        if self.pss is None: raise AttributeError("Pseudoscalar (pss) not initialized for this algebra.")
        # Return copy if value is 1 and no name, otherwise scale/name
        if value == 1 and name is None:
             # Need to decide if self.pss is mutable. Assume immutable for safety, return copy.
             return self.pss.copy() # Return a copy
        else:
             # Multiply returns a new object
             scaled_pss = self.pss * value
             # Naming symbolic MVs needs careful handling; assume name applies if symbolic
             if name is not None and scaled_pss.issymbolic:
                  # This needs a mechanism to rename symbolic components or attach name
                  log.warning(f"Naming symbolic result in pseudoscalar({name=}) not fully implemented.")
                  # scaled_pss.name = name # Placeholder - might not work
             return scaled_pss

    # --- Geometric Object Factories (Rotor, Translator, etc.) ---
    # These require specific algebra setups (e.g., PGA for translator)
    # and mathematical functions (cos, sin, sqrt). Ensure 'math' is imported.

    def rotor(self, angle: float, plane_indices: Tuple[int, int], *, name: Optional[str] = None) -> 'MultiVector':
        """Creates a rotor for rotation in a specified plane."""
        if not isinstance(plane_indices, tuple) or len(plane_indices) != 2: raise TypeError("plane_indices must be a tuple of two integers.")
        i, j = sorted(plane_indices)
        if i == j: raise ValueError("plane_indices must contain two distinct indices.")
        try:
            # R = exp(-B * angle / 2) where B is the unit bivector for the plane
            # B = blade([i, j])
            B = self.blade([i, j])
            # Check B*B to ensure it's -1 (or appropriate scalar)
            B2 = (B*B).e # Geometric product and extract scalar part
            if abs(B2 + 1.0) > 1e-9: # Check if B^2 is close to -1
                 raise ValueError(f"Cannot create rotor: Plane e{i}{j} does not square to -1 (B^2 = {B2}). Ensure Euclidean signature for these indices.")

            # Use MultiVector exp method if available, or manual calculation
            # Manual: cos(angle/2) - B * sin(angle/2)
            cos_half = math.cos(angle / 2.0)
            sin_half = math.sin(angle / 2.0)
            rotor = self.scalar(cos_half) - B * sin_half
            if name is not None and rotor.issymbolic: log.warning(f"Naming symbolic result in rotor({name=}) not fully implemented.")
            return rotor
        except (ValueError, AttributeError, TypeError) as e:
             raise ValueError(f"Failed to create rotor for plane {plane_indices}: {e}") from e


    def translator(self, direction: Sequence[float], distance: float = 1.0, *, name: Optional[str] = None) -> 'MultiVector':
        """Creates a translator in Projective Geometric Algebra (PGA). Requires r=1, start_index=0."""
        if self.r != 1 or self.start_index != 0: raise ValueError("Translator requires a PGA setup (r=1, start_index=0).")
        if len(direction) != self.d - 1: raise ValueError(f"Direction vector must have {self.d - 1} components for {self.d}D PGA.")

        # T = 1 - d/2 * n, where n is the direction vector represented as e0i, e0j, ...
        translator = self.scalar(1)
        for i, component in enumerate(direction, start=1): # Assumes non-null indices start at 1
            if component != 0:
                coeff = - (distance / 2.0) * component
                try:
                     # Blade e0i (assuming 0 is the null vector index)
                     blade_term = self.blade([0, i], value=coeff)
                     translator += blade_term
                except (ValueError, KeyError) as e:
                     raise ValueError(f"Failed to create translator blade for direction component {i}: {e}") from e
        if name is not None and translator.issymbolic: log.warning(f"Naming symbolic result in translator({name=}) not fully implemented.")
        return translator

    def reflector(self, normal: Sequence[float], *, name: Optional[str] = None) -> 'MultiVector':
        """Creates a reflector (normalized plane or vector)."""
        # Treat normal as vector for reflection. Plane representation might differ.
        if len(normal) != self.d: raise ValueError(f"Normal vector length ({len(normal)}) must match algebra dimension ({self.d}).")
        norm_vec = self.vector(values=normal, name=name)
        try:
            reflector = norm_vec.normalized()
            if name is not None and reflector.issymbolic: log.warning(f"Naming symbolic result in reflector({name=}) not fully implemented.")
            return reflector
        except ZeroDivisionError: raise ValueError("Normal vector cannot be zero (cannot normalize).")
        except TypeError as e: raise TypeError(f"Could not normalize the normal vector: {e}")


    def motor(self, rotation: Optional[Tuple[float, Tuple[int, int]]] = None, translation: Optional[Tuple[Sequence[float], float]] = None, *, name: Optional[str] = None) -> 'MultiVector':
        """Creates a motor (rotation * translation) typically in PGA."""
        if not rotation and not translation: return self.scalar(1, name=name)
        # Translation requires PGA setup
        if translation and (self.r != 1 or self.start_index != 0): raise ValueError("Translation component requires a PGA setup (r=1, start_index=0).")

        R = self.scalar(1); T = self.scalar(1)
        if rotation: angle, plane_indices = rotation; R = self.rotor(angle, plane_indices)
        if translation: direction, distance = translation; T = self.translator(direction, distance)

        # Motor M = T * R (translation followed by rotation relative to origin)
        motor = T * R
        if name is not None and motor.issymbolic: log.warning(f"Naming symbolic result in motor({name=}) not fully implemented.")
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
        if not hasattr(self, 'pss') or self.pss is None: raise AttributeError("Algebra not fully initialized or pss is None.")
        try: return self.pss.inv()
        except ZeroDivisionError as e: raise ZeroDivisionError(f"Cannot compute pss_inv: Pseudoscalar is not invertible. Original error: {e}") from e
        except Exception as e: raise RuntimeError(f"Unexpected error computing pss_inv: {e}") from e


    @property
    @lru_cache(maxsize=None) # Cache the Cayley table result
    def cayley(self) -> np.ndarray:
        """Generates the Cayley table using the compiled GP function."""
        basis_indices = np.arange(len(self))
        # Use the precompiled GP function
        vectorized_gp = np.vectorize(self._gp_func, otypes=[tuple]) # Specify otypes=tuple
        # Compute table entries
        i_indices, j_indices = np.meshgrid(basis_indices, basis_indices, indexing='ij')
        # vectorized_gp returns tuples (result_idx, sign)
        table_tuples = vectorized_gp(i_indices, j_indices)
        # If you need just the resulting indices or signs, process table_tuples
        # For now, return the array of tuples
        return table_tuples # Returns array of (result_idx, sign) tuples

    @property
    @lru_cache(maxsize=None) # Cache the frame result
    def frame(self) -> List['MultiVector']:
        """Returns a list of the grade-1 basis vectors."""
        vector_keys = self.indices_for_grades.get((1,), [])
        # Use BladeDict __getitem__ to get vector MultiVectors
        return [self.blades[self.bin2canon[key]] for key in vector_keys]

    @property
    @lru_cache(maxsize=None) # Cache the reciprocal frame result
    def reciprocal_frame(self) -> List['MultiVector']:
        """Computes the reciprocal frame vectors."""
        if self.r > 0: raise NotImplementedError("Reciprocal frame is not uniquely defined for degenerate algebras (r > 0).")
        if self.d == 0: return []
        n = self.d; basis_vectors = self.frame
        try: pss_inv_val = self.pss_inv # Get cached inverse PSS
        except (AttributeError, ZeroDivisionError, RuntimeError) as e:
            raise RuntimeError(f"Cannot compute reciprocal frame: Failed to get PSS inverse. {e}") from e

        reciprocal_vectors = []
        all_vector_indices = list(range(n))
        # Compute using formula: e^i = (-1)^(i*(n-1)) * (e_1 ^ ... ^ e_{i-1} ^ e_{i+1} ^ ... ^ e_n) * I^{-1}
        # Equivalent to (-1)^i * wedge_product_without_i * pss_inv
        for i in range(n):
            wedge_indices = [j for j in all_vector_indices if j != i]
            if not wedge_indices: # Should only happen if d=1
                 wedge_product = self.scalar(1)
            else:
                 # Compute wedge product of remaining vectors
                 wedge_product = reduce(operator.xor, [basis_vectors[j] for j in wedge_indices])
            sign = (-1)**i
            reciprocal_vec = sign * wedge_product * pss_inv_val
            reciprocal_vectors.append(reciprocal_vec)
        return reciprocal_vectors

    def register(self, _func: Optional[Callable] = None, *, symbolic: bool = False, **options) -> Callable:
        """Decorator to register a custom function with the algebra."""
        # (Implementation unchanged)
        def decorator_register(func: Callable) -> Callable:
            @wraps(func)
            def wrapper(*args, **kwargs):
                log.info(f"Calling registered function '{func.__name__}' (symbolic={symbolic}, options={options})")
                # Potentially add logic here to use codegen if symbolic=True?
                return func(*args, **kwargs)
            log.info(f"Registering function '{func.__name__}' with algebra {self}. Options: symbolic={symbolic}, {options}")
            self.registry[func.__name__] = wrapper
            # Allow access via attribute, e.g., alg.my_custom_func(...)
            setattr(self, func.__name__, wrapper)
            return wrapper
        if _func is None: return decorator_register
        else: return decorator_register(_func)

    def graph(self, *args, **kwargs) -> Any:
        """Creates a GraphWidget for visualization."""
        # Import locally to avoid hard dependency if graph extras aren't installed
        try:
            from kingdon.graph import GraphWidget
            from kingdon.multivector import MultiVector # Import locally too
        except ImportError: log.error("GraphWidget or MultiVector could not be imported. Please ensure 'kingdon[graph]' extras are installed."); raise
        # Filter args to pass only MultiVectors to GraphWidget constructor if needed
        # Or pass all args and let GraphWidget handle them
        return GraphWidget(algebra=self, raw_subjects=list(args), options=kwargs)

    def __len__(self) -> int: return 2**self.d

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, Algebra): return NotImplemented
        sig_equal = np.array_equal(self.signature, other.signature) if self.signature is not None and other.signature is not None else self.signature == other.signature
        basis_names_equal = self.basis_names == other.basis_names
        # Include codegen_symbolcls comparison? Maybe too strict.
        return (self.p == other.p and self.q == other.q and self.r == other.r and
                self.start_index == other.start_index and sig_equal and
                basis_names_equal and self.graded == other.graded)

    def __hash__(self) -> int:
        sig_tuple = tuple(self.signature.tolist()) if self.signature is not None else None
        basis_tuple = tuple(self.basis_names) if self.basis_names else None
        # Include relevant state in hash
        return hash((self.p, self.q, self.r, sig_tuple, self.start_index, basis_tuple, self.graded))

    def grade(self, mv: 'MultiVector', *grades: int) -> 'MultiVector':
         """Extracts specific grades from a multivector."""
         # Delegate to multivector's grade method
         return mv.grade(*grades)

    def dual(self, mv: 'MultiVector', kind: str = 'auto') -> 'MultiVector':
         """Computes the dual of a multivector."""
         # Delegate to multivector's dual method
         return mv.dual(kind=kind)

    def undual(self, mv: 'MultiVector', kind: str = 'auto') -> 'MultiVector':
         """Computes the undual of a multivector."""
         # Delegate to multivector's undual method
         return mv.undual(kind=kind)

    @classmethod
    def from_signature(cls, signature: List[int], start_index: Optional[int] = None) -> 'Algebra':
         """Creates an Algebra instance from an explicit signature list."""
         return cls(signature=signature, start_index=start_index)

    @classmethod
    def from_name(cls, name: str) -> 'Algebra':
        """Creates an Algebra instance from a common name (e.g., 'PGA3D')."""
        # (Implementation unchanged)
        name_upper = name.upper().strip()
        # Define common algebras
        algebras = {
             "CGA2D": {'p': 3, 'q': 1, 'r': 0, 'start_index': 1},
             "CGA3D": {'p': 4, 'q': 1, 'r': 0, 'start_index': 1},
             "PGA2D": {'p': 2, 'q': 0, 'r': 1, 'start_index': 0},
             "PGA3D": {'p': 3, 'q': 0, 'r': 1, 'start_index': 0},
             "3DPGA": {'p': 3, 'q': 0, 'r': 1, 'start_index': 0},
             "VGA2D": {'p': 2, 'q': 0, 'r': 0, 'start_index': 1},
             "VGA3D": {'p': 3, 'q': 0, 'r': 0, 'start_index': 1},
             "STA":   {'p': 1, 'q': 3, 'r': 0, 'start_index': 0}, # Spacetime Algebra (t, x, y, z)
             "APS":   {'p': 3, 'q': 0, 'r': 0, 'start_index': 1}  # Alias for VGA3D
        }
        if name_upper in algebras: return cls(**algebras[name_upper])
        else:
             # Support Cl(p,q,r) or VGA(p,q,r) syntax
             match_cl = re.match(r'^(?:CL|VGA)\((\d+)(?:,(\d+))?(?:,(\d+))?\)$', name_upper)
             if match_cl:
                 p = int(match_cl.group(1))
                 q = int(match_cl.group(2) or 0)
                 r = int(match_cl.group(3) or 0)
                 # Default start_index for Cl(p,q,r) - use 1 unless r=1
                 start_idx = 0 if r == 1 else 1
                 return cls(p=p, q=q, r=r, start_index=start_idx)
             raise ValueError(f"Unknown or invalid algebra name: {name}")

# Convenience factory functions calling Algebra.from_name (Unchanged)
def PGA2D() -> Algebra: return Algebra.from_name("PGA2D")
def PGA3D() -> Algebra: return Algebra.from_name("PGA3D")
def CGA2D() -> Algebra: return Algebra.from_name("CGA2D")
def CGA3D() -> Algebra: return Algebra.from_name("CGA3D")
def VGA3D() -> Algebra: return Algebra.from_name("VGA3D")
def STA() -> Algebra: return Algebra.from_name("STA")
APS = VGA3D # Alias