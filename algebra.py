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


from .operator_dict import OperatorDict, UnaryOperatorDict, BladeDict, Registry, AlgebraError
from .multivector   import MultiVector
from .codegen       import CodegenOutput

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

    # You need to ADD this method to the Algebra class in the new file.

    def _initialize_blade_dictionary(self, BladeDict_cls, MultiVector_cls):
        """
        Initializes the dictionary of basis blades and the pss attribute.
        """
        # Blades are not precomputed for algebras larger than 6D.
        # Note: The original BladeDict from Originalalgebra.py is what is needed here.
        self.blades = BladeDict_cls(algebra=self, lazy=self.d > 6)

        # The pss is the blade with the highest grade.
        pss_key = (1 << self.d) - 1
        self.pss = self.blades[self.bin2canon[pss_key]]
    
    
    # ADD THIS METHOD to the Algebra class in algebra.py

    def _canonical_key_sort_func(self) -> Callable[[int], Tuple[int, float]]:
        """
        Returns a callable that can be used as a sort key to sort basis blade
        keys in canonical order (by grade, then by value).
        """
        # The key function uses the pre-computed map for efficiency.
        # We return the .get method of the dict, which is a callable.
        return self._key_to_grade_index_map.get
    
    
    # In algebra.py, inside the Algebra class

    def _parse_composite_custom_key(self, basis_blade: str) -> bool:
        """
        Placeholder for parsing composite custom basis names.
        Currently returns False, indicating it did not handle the key.
        """
        # This is a placeholder. A full implementation would check if `basis_blade`
        # can be formed by concatenating names from `self.basis_names`.
        return False
    
    # ADD THIS METHOD to the Algebra class in algebra.py

    def _compile_gp_func(self) -> Callable[[int, int], Tuple[int, int]]:
        """
        Creates and returns a JIT-compilable or simple Python function for the
        geometric product of two basis blade keys.
        
        This function is a performance-critical part for operations like
        building the Cayley table.
        """
        # We can use a simple lambda here because the self.signs dictionary
        # is already pre-computed and available in the closure.
        # It returns a tuple of (resulting_key, sign).
        gp_func = lambda i, j: (i ^ j, self.signs.get((i, j), 0))

        # If a JIT wrapper like Numba is available on the algebra, wrap it.
        # This provides a significant performance boost.
        if self.wrapper and callable(self.wrapper):
            try:
                # Attempt to wrap the function. This might fail if the wrapper
                # (e.g., Numba) doesn't support dict lookups, but for basic
                # operations, this pattern is common. If it fails, we fall
                # back to the plain Python version.
                return self.wrapper(gp_func)
            except Exception as e:
                log.warning(f"Could not apply JIT wrapper to gp_func, falling back to pure Python. Error: {e}")
                return gp_func
        else:
            # Return the plain Python function if no wrapper is specified.
            return gp_func
    
    # ADD THIS METHOD to the Algebra class in algebra.py

    def _create_registry(self):
        """
        Initializes the operator dictionaries (e.g., gp, sw) and populates
        the main registry.
        """
        from . import codegen  # Local import to get codegen functions

        self.registry = {}
        for op_name, meta in self._operator_metadata.items():
            # Get the codegen function from the codegen module by its string name
            codegen_func = getattr(codegen, meta['codegen'], None)
            if codegen_func is None:
                log.warning(f"Codegen function '{meta['codegen']}' for operator '{op_name}' not found. Skipping.")
                continue

            # Get the OperatorDict class (e.g., OperatorDict, UnaryOperatorDict)
            op_cls = meta['cls']

            # Create the instance
            op_instance = op_cls(name=op_name, codegen=codegen_func, algebra=self)
            
            # Set the operator on the algebra instance (e.g., self.gp = ...)
            setattr(self, op_name, op_instance)
            # Add it to the main registry dictionary
            self.registry[op_name] = op_instance
    
    
    # Add _grade_index_map caching to __post_init__ for efficiency
    def __post_init__(self) -> None:
        # Imports moved here to avoid potential circularity at module level
        # Need to handle potential import errors if modules aren't found
        try:
            from .multivector import MultiVector
            from .operator_dict import BladeDict
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
        
        # Explicitly handle e0 mapping if starting at 0, which might not be a basis vector
        if self.start_index == 0 and "e0" not in self.canon2bin:
            # Find which basis vector corresponds to index 0
            zero_idx_vector_key = 1 << (0 - self.start_index)
            if self.signature[0] == 0: # If e0 is null, it's special
                self.canon2bin["e0"] = zero_idx_vector_key
            # Else, if it's not null, 'e0' can just be an alias for the first vector
            # However, the universal alias in _blade2canon handles 'e0' -> scalar, which is often desired.
            # This part of the logic is subtle. The key is that `_blade2canon` handles the *input*,
            # and this part handles the *canonical representation*.

    def _precompute_grade_indices(self) -> None:
        """
        Fills self.indices_for_grades with all basis indices for each grade.
        """
        self.indices_for_grades = defaultdict(list)
        all_keys_iterable = range(len(self))
        key_func = lambda k: bin(k).count('1')
        all_keys_sorted = sorted(all_keys_iterable, key=key_func)
        for grade, group in groupby(all_keys_sorted, key=key_func):
            self.indices_for_grades[(grade,)] = list(group)

    def _compute_multiplication_signs(self) -> None:
        """
        Precomputes the sign of the geometric product for all pairs of basis blades.
        The result is stored in self.signs.
        """
        self.signs = {};
        algebra_len = len(self)
        for i in range(algebra_len):
            for j in range(algebra_len):
                if i == 0 or j == 0:
                    self.signs[(i, j)] = 1
                    continue
                
                # Decompose basis blades into their vector indices
                i_indices = {k for k in range(self.d) if (i >> k) & 1}
                j_indices = {k for k in range(self.d) if (j >> k) & 1}
                
                common_indices = i_indices.intersection(j_indices)
                i_only_indices = i_indices.difference(common_indices)
                
                # Count swaps required to move all elements of i_only_indices past j_indices
                swaps = sum(1 for i_idx in i_only_indices for j_idx in j_indices if i_idx > j_idx)
                
                sign = (-1)**swaps
                
                # Multiply by the squares of common basis vectors
                for k in common_indices:
                    sig_val = self.signature[k]
                    if sig_val == 0:
                        sign = 0; break
                    sign *= sig_val
                
                if sign != 0:
                    self.signs[(i, j)] = int(sign)

    # This method replaces the old `_blade2canon` completely.
    def _blade2canon(self, basis_blade: str) -> Tuple[str, int]:
        """
        Converts a user-provided blade string into its canonical name and sign swap count.
        Handles 'e0' as a universal alias for the scalar 'e'.
        """
        if not isinstance(basis_blade, str):
            raise TypeError(f"Blade names must be strings, got {type(basis_blade).__name__}.")

        # Handle universal scalar aliases first.
        if basis_blade == "e" or basis_blade == "e0":
            return "e", 0

        # Check for direct match in canonical names (e.g., 'e1', 'e12', or custom names)
        if basis_blade in self.canon2bin:
            return basis_blade, 0

        # If not a direct match, parse as `e<indices>` format
        match_default = re.match(r'^e([0-9]+)$', basis_blade)
        if not match_default:
            # Try to parse as a composite of custom names if available
            if self.basis_names and self._parse_composite_custom_key(basis_blade):
                return basis_blade, 0 # Let the creation logic handle composites
            raise ValueError(f"Invalid basis blade format: '{basis_blade}'. Not a canonical name or 'e<indices>'.")

        indices_str = match_default.group(1)
        indices = [int(c) for c in indices_str]

        if len(set(indices)) != len(indices):
            raise ValueError(f"Duplicate indices in blade name: {basis_blade}")

        min_idx, max_idx = self.start_index, self.start_index + self.d - 1
        if not all(min_idx <= i <= max_idx for i in indices):
            raise ValueError(f"Indices in '{basis_blade}' must be between {min_idx} and {max_idx}")

        sorted_indices = sorted(indices)
        canonical_name = f"e{''.join(map(str, sorted_indices))}"
        
        swaps = 0
        current_indices = list(indices)
        n = len(current_indices)
        for i in range(n):
            for j in range(n - 1 - i):
                if current_indices[j] > current_indices[j+1]:
                    current_indices[j], current_indices[j+1] = current_indices[j+1], current_indices[j]
                    swaps += 1
        return canonical_name, swaps

    # This method is now DELETED. It is replaced by the logic within the new `multivector` method.
    # def _process_creation_args(...):
    #     ...

    # This is the new, definitive `multivector` factory method.
    # It replaces the old `multivector` method and the `_process_creation_args` logic.
    def multivector(self, values: Optional[Any] = None, *,
                    keys: Optional[Sequence[Union[int, str]]] = None,
                    grades: Optional[Sequence[int]] = None,
                    name: Optional[str] = None,
                    **blade_kwargs: Any) -> 'MultiVector':
        """
        Primary factory method to create a multivector in this algebra.

        Order of precedence: `name` (symbolic), `keys`, `grades`, `blade_kwargs`, `values`.
        """
        from .multivector import MultiVector # Local import

        # --- Case 1: Symbolic Creation ---
        if name is not None:
            if blade_kwargs or values is not None:
                raise ValueError("Cannot specify 'values' or blade keywords with 'name' for symbolic creation.")
            
            sym_cls = self.codegen_symbolcls
            keys_tuple: Tuple[int, ...]

            if keys is not None:
                key_set = {self.canon2bin[self._blade2canon(k)[0]] if isinstance(k, str) else k for k in keys}
                keys_tuple = tuple(sorted(list(key_set), key=self._canonical_key_sort_func()))
            elif grades is not None:
                key_set = {k for g in grades for k in self.indices_for_grades.get((g,), [])}
                keys_tuple = tuple(sorted(list(key_set), key=self._canonical_key_sort_func()))
            else:
                # Full symbolic multivector if no keys/grades specified
                keys_tuple = tuple(range(len(self)))

            final_values = [sym_cls(f"{name}_{self.bin2canon.get(k, k)}") for k in keys_tuple]
            return MultiVector.fromkeysvalues(self, keys=keys_tuple, values=final_values)

        # --- Case 2: Creation from Blade Keywords ---
        if blade_kwargs:
            if values is not None or keys is not None or grades is not None:
                raise ValueError("Cannot mix blade keywords with 'values', 'keys', or 'grades'.")
            
            key_set = {self.canon2bin[self._blade2canon(k)[0]] for k in blade_kwargs.keys()}
            keys_tuple = tuple(sorted(list(key_set), key=self._canonical_key_sort_func()))
            # Map canonical names back to binary keys to retrieve values correctly
            canon_map = {self.canon2bin[self._blade2canon(k)[0]]: v for k, v in blade_kwargs.items()}
            final_values = [canon_map.get(k, 0) for k in keys_tuple]
            return MultiVector.fromkeysvalues(self, keys=keys_tuple, values=final_values)

        # --- Case 3: Creation from explicit keys and values ---
        if keys is not None:
            key_set = {self.canon2bin[self._blade2canon(k)[0]] if isinstance(k, str) else k for k in keys}
            keys_tuple = tuple(sorted(list(key_set), key=self._canonical_key_sort_func()))
            return MultiVector.fromkeysvalues(self, keys=keys_tuple, values=values)
            
        # --- Case 4: Creation from grades and values ---
        if grades is not None:
            key_set = {k for g in grades for k in self.indices_for_grades.get((g,), [])}
            keys_tuple = tuple(sorted(list(key_set), key=self._canonical_key_sort_func()))
            if values is None:
                values = [1] * len(keys_tuple)
            return MultiVector.fromkeysvalues(self, keys=keys_tuple, values=values)

        # --- Case 5: Creation from a dictionary of blade names to values ---
        if isinstance(values, dict):
            # Handle empty dictionary case - create truly empty multivector
            if not values:
                return MultiVector.fromkeysvalues(self, keys=(), values=[])
            # If all keys are ints, treat as blade-indexed values (the main GA case!)
            if all(isinstance(k, int) for k in values.keys()):
                keys_tuple = tuple(sorted(values.keys()))
                final_values = [values[k] for k in keys_tuple]
                return MultiVector.fromkeysvalues(self, keys=keys_tuple, values=final_values)
            # If keys are strings, treat as blade keywords
            elif all(isinstance(k, str) for k in values.keys()):
                return self.multivector(**values)
            else:
                raise TypeError("Dictionary keys for 'values' must be all ints (for blade indices) or all strings (for blade keywords).")

        # --- Case 6: Fallback for scalar, full list/array, or empty ---
        if values is None: values = 0
        if isinstance(values, (int, float, complex, sympy.Expr)):
            return MultiVector.fromkeysvalues(self, keys=(0,), values=[values])
        if isinstance(values, (list, tuple, np.ndarray)):
            if len(values) != len(self):
                raise ValueError(f"Ambiguous list `values` of length {len(values)}. Must match algebra length {len(self)} or provide `keys`/`grades`.")
            final_keys_values = [(i, v) for i, v in enumerate(values) if v != 0]
            if not final_keys_values: return self.scalar(0)
            final_keys, final_values = zip(*final_keys_values)
            return MultiVector.fromkeysvalues(self, keys=final_keys, values=list(final_values))
        
        raise TypeError(f"Could not create MultiVector from supplied arguments: {values}")

    # ========================================================================
    # Convenience Factory Methods (Now calling the new self.multivector)
    # ========================================================================

    def scalar(self, value: Any = 0, *, name: Optional[str] = None) -> 'MultiVector':
        """Creates a scalar (grade 0) multivector."""
        return self.multivector(values=value, name=name, grades=(0,))

    def vector(self, values: Optional[Any] = None, *, name: Optional[str] = None) -> 'MultiVector':
        """Creates a vector (grade 1) multivector."""
        return self.multivector(values=values, name=name, grades=(1,))

    def bivector(self, values: Optional[Any] = None, *, name: Optional[str] = None) -> 'MultiVector':
        """Creates a bivector (grade 2) multivector."""
        return self.multivector(values=values, name=name, grades=(2,))

    def trivector(self, values: Optional[Any] = None, *, name: Optional[str] = None) -> 'MultiVector':
        """Creates a trivector (grade 3) multivector."""
        return self.multivector(values=values, name=name, grades=(3,))

    def purevector(self, grade: int, values: Optional[Any] = None, *, name: Optional[str] = None) -> 'MultiVector':
        """Creates a pure grade multivector."""
        if not isinstance(grade, int) or not (0 <= grade <= self.d):
            raise ValueError(f"Grade must be an integer between 0 and {self.d}")
        return self.multivector(values=values, grades=(grade,), name=name)

    def evenmv(self, values: Optional[Any] = None, *, name: Optional[str] = None) -> 'MultiVector':
        """Creates a multivector containing only even grades."""
        even_grades = tuple(g for g in range(0, self.d + 1) if g % 2 == 0)
        return self.multivector(values=values, grades=even_grades, name=name)

    def oddmv(self, values: Optional[Any] = None, *, name: Optional[str] = None) -> 'MultiVector':
        """Creates a multivector containing only odd grades."""
        odd_grades = tuple(g for g in range(1, self.d + 1) if g % 2 != 0)
        return self.multivector(values=values, grades=odd_grades, name=name)

    def blade(self, spec: Union[str, Sequence[int]], value: Any = 1, *, name: Optional[str] = None) -> 'MultiVector':
        """Creates a multivector representing a single basis blade."""
        blade_values = {}
        sign = 1
        if isinstance(spec, (list, tuple, set)):
            indices = list(spec)
            if not indices: key_str = "e"
            else:
                min_idx, max_idx = self.start_index, self.start_index + self.d - 1
                if not all(isinstance(i, int) and min_idx <= i <= max_idx for i in indices):
                    raise ValueError(f"Indices {indices} must be integers between {min_idx} and {max_idx}")
                if len(set(indices)) != len(indices): raise ValueError(f"Indices in {indices} must be unique.")
                sorted_indices = sorted(indices)
                key_str = f"e{''.join(map(str, sorted_indices))}"
                current_indices = list(indices); swaps = 0; n = len(current_indices)
                for i in range(n):
                    for j in range(n - 1 - i):
                        if current_indices[j] > current_indices[j+1]: current_indices[j], current_indices[j+1] = current_indices[j+1], current_indices[j]; swaps += 1
                sign = (-1)**swaps
        elif isinstance(spec, str):
            canonical_name, swaps = self._blade2canon(spec)
            key_str = canonical_name
            sign = (-1)**swaps
        else: raise TypeError("Blade specification must be a string or an iterable of indices.")

        final_value = sign * value
        blade_values[key_str] = final_value
        return self.multivector(name=name, **blade_values)

    def pseudoscalar(self, value: Any = 1, *, name: Optional[str] = None) -> 'MultiVector':
        """Creates the pseudoscalar multivector, optionally scaled."""
        if self.pss is None: raise AttributeError("Pseudoscalar (pss) not initialized for this algebra.")
        if value == 1 and name is None:
            return self.pss.copy()
        else:
            scaled_pss = self.pss * value
            if name is not None:
                # This is tricky for symbolic results. The best we can do is create a new symbolic
                # object if the result is symbolic. For now, this is a known limitation.
                log.warning(f"Naming the result of a scaled pseudoscalar is not fully supported.")
            return scaled_pss

    def rotor(self, angle: float, plane_indices: Tuple[int, int], *, name: Optional[str] = None) -> 'MultiVector':
        """Creates a rotor for rotation in a specified plane."""
        if not isinstance(plane_indices, tuple) or len(plane_indices) != 2:
            raise TypeError("plane_indices must be a tuple of two integers.")
        B = self.blade(plane_indices)
        B2_mv = B*B
        if not (B2_mv.grades == (0,) and B2_mv.e < 0):
            raise ValueError(f"Plane {plane_indices} does not square to a negative scalar. B*B = {B2_mv}")
        
        cos_half = math.cos(angle / 2.0)
        sin_half = math.sin(angle / 2.0)
        
        rotor = self.scalar(cos_half) + (B / math.sqrt(-B2_mv.e)) * sin_half
        return rotor # Naming is complex for composite objects, deferring.

    def translator(self, direction: Sequence[float], distance: float = 1.0, *, name: Optional[str] = None) -> 'MultiVector':
        """Creates a translator in Projective Geometric Algebra (PGA)."""
        if self.r != 1 or self.start_index != 0:
            raise ValueError("Translator requires a PGA setup (r=1, start_index=0).")
        if len(direction) != self.d - 1:
            raise ValueError(f"Direction vector must have {self.d - 1} components for {self.d}D PGA.")

        translator = self.scalar(1.0)
        for i, component in enumerate(direction, start=1):
            if component != 0:
                coeff = - (distance / 2.0) * component
                translator += self.blade([0, i], value=coeff)
        return translator

    def reflector(self, normal: Sequence[float], *, name: Optional[str] = None) -> 'MultiVector':
        """Creates a reflector (normalized plane or vector)."""
        if len(normal) > self.d:
            raise ValueError(f"Normal vector length ({len(normal)}) cannot exceed algebra dimension ({self.d}).")
        
        # Pad with zeros if a lower-dimensional vector is given
        padded_normal = list(normal) + [0] * (self.d - len(normal))
        
        norm_vec = self.vector(values=padded_normal, name=name)
        if norm_vec.norm() == 0:
            raise ValueError("Normal vector cannot be zero.")
        return norm_vec.normalized()

    def motor(self, rotation: Optional[Tuple[float, Tuple[int, int]]] = None, translation: Optional[Tuple[Sequence[float], float]] = None, *, name: Optional[str] = None) -> 'MultiVector':
        """Creates a motor (rotation * translation) typically in PGA."""
        R = self.scalar(1)
        if rotation:
            angle, plane_indices = rotation
            R = self.rotor(angle, plane_indices)
        
        T = self.scalar(1)
        if translation:
            direction, distance = translation
            T = self.translator(direction, distance)
            
        motor = T * R
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
            from .graph import GraphWidget
            from .multivector import MultiVector # Import locally too
        except ImportError: log.error("GraphWidget or MultiVector could not be imported. Please ensure 'kingdon[graph]' extras are installed."); raise
        # Filter args to pass only MultiVectors to GraphWidget constructor if needed
        # Or pass all args and let GraphWidget handle them
        return GraphWidget(algebra=self, raw_subjects=list(args), options=kwargs)

    def __len__(self) -> int: return 2**self.d

    

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

    def gp(self, x: 'MultiVector', y: 'MultiVector') -> 'MultiVector':
        """
        Compute the geometric product of two multivectors.
        
        Args:
            x: First multivector
            y: Second multivector
            
        Returns:
            Geometric product x * y
        """
        from .multivector import MultiVector
        
        if x.algebra != self or y.algebra != self:
            raise ValueError("Multivectors must belong to this algebra")
            
        # Result dictionary to accumulate terms
        result_terms = {}
        
        # Iterate over all combinations of basis blades
        for i, (kx, vx) in enumerate(zip(x._keys, x._values)):
            for j, (ky, vy) in enumerate(zip(y._keys, y._values)):
                # Get the sign from the Cayley table
                sign = self.signs[kx, ky] if hasattr(self, 'signs') else 1
                
                if sign != 0:  # Non-zero product
                    # XOR the keys to get output key (geometric product rule)
                    key_out = kx ^ ky
                    
                    # Compute the coefficient
                    if isinstance(vx, (int, float, complex)) and isinstance(vy, (int, float, complex)):
                        coeff = sign * vx * vy
                    else:
                        # Handle symbolic case
                        import sympy
                        coeff = sign * sympy.Mul(vx, vy, evaluate=False) if sign > 0 else -sympy.Mul(vx, vy, evaluate=False)
                    
                    # Accumulate terms
                    if key_out in result_terms:
                        result_terms[key_out] += coeff
                    else:
                        result_terms[key_out] = coeff
        
        # Filter out zero terms and create result
        if not result_terms:
            return MultiVector.fromkeysvalues(self, keys=(), values=[])
            
        # Sort keys and extract values
        sorted_keys = tuple(sorted(result_terms.keys()))
        sorted_values = [result_terms[k] for k in sorted_keys]
        
        return MultiVector.fromkeysvalues(self, keys=sorted_keys, values=sorted_values)

    def inner_product(self, x: 'MultiVector', y: 'MultiVector') -> 'MultiVector':
        """
        Compute the inner product of two multivectors.
        
        Args:
            x: First multivector
            y: Second multivector
            
        Returns:
            Inner product x | y
        """
        return self._contraction_product(x, y, diff_func=abs)
    
    def left_contraction(self, x: 'MultiVector', y: 'MultiVector') -> 'MultiVector':
        """
        Compute the left contraction of two multivectors.
        
        Args:
            x: First multivector
            y: Second multivector
            
        Returns:
            Left contraction x << y
        """
        return self._contraction_product(x, y, diff_func=lambda k: -k)
    
    def right_contraction(self, x: 'MultiVector', y: 'MultiVector') -> 'MultiVector':
        """
        Compute the right contraction of two multivectors.
        
        Args:
            x: First multivector
            y: Second multivector
            
        Returns:
            Right contraction x >> y
        """
        return self._contraction_product(x, y, diff_func=lambda k: k)
    
    def _contraction_product(self, x: 'MultiVector', y: 'MultiVector', diff_func) -> 'MultiVector':
        """
        Helper method for computing contraction-based products.
        
        Args:
            x: First multivector
            y: Second multivector
            diff_func: Function to apply to key difference (abs for inner, -k for left, k for right)
            
        Returns:
            Contraction product result
        """
        from .multivector import MultiVector
        
        if x.algebra != self or y.algebra != self:
            raise ValueError("Multivectors must belong to this algebra")
            
        # Result dictionary to accumulate terms
        result_terms = {}
        
        # Iterate over all combinations of basis blades
        for kx, vx in zip(x._keys, x._values):
            for ky, vy in zip(y._keys, y._values):
                # Get the sign from the Cayley table
                sign = self.signs[kx, ky] if hasattr(self, 'signs') else 1
                
                if sign != 0:  # Non-zero product
                    # Apply contraction filter: k_out == diff_func(kx - ky)
                    key_diff = kx - ky
                    expected_key_out = diff_func(key_diff)
                    actual_key_out = kx ^ ky
                    
                    # Only include terms that satisfy the contraction condition
                    if actual_key_out == expected_key_out:
                        # Compute the coefficient
                        if isinstance(vx, (int, float, complex)) and isinstance(vy, (int, float, complex)):
                            coeff = sign * vx * vy
                        else:
                            # Handle symbolic case
                            import sympy
                            coeff = sign * sympy.Mul(vx, vy, evaluate=False) if sign > 0 else -sympy.Mul(vx, vy, evaluate=False)
                        
                        # Accumulate terms
                        if actual_key_out in result_terms:
                            result_terms[actual_key_out] += coeff
                        else:
                            result_terms[actual_key_out] = coeff
        
        # Filter out zero terms and create result
        if not result_terms:
            return MultiVector.fromkeysvalues(self, keys=(), values=[])
            
        # Sort keys and extract values
        sorted_keys = tuple(sorted(result_terms.keys()))
        sorted_values = [result_terms[k] for k in sorted_keys]
        
        return MultiVector.fromkeysvalues(self, keys=sorted_keys, values=sorted_values)

    # Aliases for compatibility with existing code
    def lc(self, x: 'MultiVector', y: 'MultiVector') -> 'MultiVector':
        """Alias for left_contraction."""
        return self.left_contraction(x, y)
    
    def rc(self, x: 'MultiVector', y: 'MultiVector') -> 'MultiVector':
        """Alias for right_contraction."""
        return self.right_contraction(x, y)
    
    def ip(self, x: 'MultiVector', y: 'MultiVector') -> 'MultiVector':
        """Alias for inner_product."""
        return self.inner_product(x, y)
    
    def outer_product(self, x: 'MultiVector', y: 'MultiVector') -> 'MultiVector':
        """
        Compute the outer (exterior) product of two multivectors.
        
        Args:
            x: First multivector
            y: Second multivector
            
        Returns:
            Outer product x ^ y
        """
        from .multivector import MultiVector
        
        if x.algebra != self or y.algebra != self:
            raise ValueError("Multivectors must belong to this algebra")
            
        # Result dictionary to accumulate terms
        result_terms = {}
        
        # Iterate over all combinations of basis blades
        for kx, vx in zip(x._keys, x._values):
            for ky, vy in zip(y._keys, y._values):
                # Get the sign from the Cayley table
                sign = self.signs[kx, ky] if hasattr(self, 'signs') else 1
                
                if sign != 0:  # Non-zero product
                    # Outer product filter: k_out == kx + ky (grade addition)
                    key_out = kx ^ ky
                    grade_x = bin(kx).count('1')
                    grade_y = bin(ky).count('1')
                    grade_out = bin(key_out).count('1')
                    
                    # Only include terms where grades add up (outer product condition)
                    if grade_out == grade_x + grade_y:
                        # Compute the coefficient
                        if isinstance(vx, (int, float, complex)) and isinstance(vy, (int, float, complex)):
                            coeff = sign * vx * vy
                        else:
                            # Handle symbolic case
                            import sympy
                            coeff = sign * sympy.Mul(vx, vy, evaluate=False) if sign > 0 else -sympy.Mul(vx, vy, evaluate=False)
                        
                        # Accumulate terms
                        if key_out in result_terms:
                            result_terms[key_out] += coeff
                        else:
                            result_terms[key_out] = coeff
        
        # Filter out zero terms and create result
        if not result_terms:
            return MultiVector.fromkeysvalues(self, keys=(), values=[])
            
        # Sort keys and extract values
        sorted_keys = tuple(sorted(result_terms.keys()))
        sorted_values = [result_terms[k] for k in sorted_keys]
        
        return MultiVector.fromkeysvalues(self, keys=sorted_keys, values=sorted_values)
    
    def op(self, x: 'MultiVector', y: 'MultiVector') -> 'MultiVector':
        """Alias for outer_product."""
        return self.outer_product(x, y)
    
    def sandwich_product(self, x: 'MultiVector', y: 'MultiVector') -> 'MultiVector':
        """
        Compute the sandwich product of two multivectors: x * y * ~x.
        
        Args:
            x: First multivector (the "sandwiching" multivector)
            y: Second multivector (the multivector being sandwiched)
            
        Returns:
            Sandwich product x * y * ~x
        """
        # Sandwich product is x * y * reverse(x)
        x_rev = x.reverse()
        return self.gp(self.gp(x, y), x_rev)
    
    def sw(self, x: 'MultiVector', y: 'MultiVector') -> 'MultiVector':
        """Alias for sandwich_product."""
        return self.sandwich_product(x, y)
    
    def regressive_product(self, x: 'MultiVector', y: 'MultiVector') -> 'MultiVector':
        """
        Compute the regressive product (meet) of two multivectors.
        
        Args:
            x: First multivector
            y: Second multivector
            
        Returns:
            Regressive product x ∨ y
        """
        from .multivector import MultiVector
        
        if x.algebra != self or y.algebra != self:
            raise ValueError("Multivectors must belong to this algebra")
        
        # Regressive product is dual of outer product of duals
        # x ∨ y = *((*x) ∧ (*y))
        try:
            x_dual = x.dual()
            y_dual = y.dual()
            outer_result = self.outer_product(x_dual, y_dual)
            return outer_result.dual()
        except (AttributeError, NotImplementedError):
            # If dual is not implemented, use direct computation
            return self._regressive_product_direct(x, y)
    
    def _regressive_product_direct(self, x: 'MultiVector', y: 'MultiVector') -> 'MultiVector':
        """
        Direct computation of regressive product without using duals.
        Based on the original kingdon implementation.
        """
        from .multivector import MultiVector
        
        # Get pseudoscalar key (all bits set)
        key_pss = len(self) - 1
        
        # Result dictionary to accumulate terms
        result_terms = {}
        
        # Iterate over all combinations of basis blades
        for kx, vx in zip(x._keys, x._values):
            for ky, vy in zip(y._keys, y._values):
                # Regressive product conditions
                key_out = key_pss - (kx ^ ky)
                
                # Filter condition: key_pss == kx + ky - key_out
                if key_pss == kx + ky - key_out:
                    # Complex sign calculation for regressive product
                    if hasattr(self, 'signs'):
                        sign = (
                            self.signs[kx, key_pss - kx] *
                            self.signs[ky, key_pss - ky] *
                            self.signs[key_pss - kx, key_pss - ky] *
                            self.signs[key_pss - (kx ^ ky), kx ^ ky]
                        )
                    else:
                        sign = 1
                    
                    if sign != 0:
                        # Compute the coefficient
                        if isinstance(vx, (int, float, complex)) and isinstance(vy, (int, float, complex)):
                            coeff = sign * vx * vy
                        else:
                            # Handle symbolic case
                            import sympy
                            coeff = sign * sympy.Mul(vx, vy, evaluate=False) if sign > 0 else -sympy.Mul(vx, vy, evaluate=False)
                        
                        # Accumulate terms
                        if key_out in result_terms:
                            result_terms[key_out] += coeff
                        else:
                            result_terms[key_out] = coeff
        
        # Filter out zero terms and create result
        if not result_terms:
            return MultiVector.fromkeysvalues(self, keys=(), values=[])
            
        # Sort keys and extract values
        sorted_keys = tuple(sorted(result_terms.keys()))
        sorted_values = [result_terms[k] for k in sorted_keys]
        
        return MultiVector.fromkeysvalues(self, keys=sorted_keys, values=sorted_values)
    
    def rp(self, x: 'MultiVector', y: 'MultiVector') -> 'MultiVector':
        """Alias for regressive_product."""
        return self.regressive_product(x, y)
    
    def commutator_product(self, x: 'MultiVector', y: 'MultiVector') -> 'MultiVector':
        """
        Compute the commutator product of two multivectors: 0.5 * (x*y - y*x).
        
        Args:
            x: First multivector
            y: Second multivector
            
        Returns:
            Commutator product [x, y] = 0.5 * (x*y - y*x)
        """
        xy = self.gp(x, y)
        yx = self.gp(y, x)
        diff = xy - yx
        return diff * 0.5
    
    def cp(self, x: 'MultiVector', y: 'MultiVector') -> 'MultiVector':
        """Alias for commutator_product."""
        return self.commutator_product(x, y)
    
    def anticommutator_product(self, x: 'MultiVector', y: 'MultiVector') -> 'MultiVector':
        """
        Compute the anti-commutator product of two multivectors: 0.5 * (x*y + y*x).
        
        Args:
            x: First multivector
            y: Second multivector
            
        Returns:
            Anti-commutator product {x, y} = 0.5 * (x*y + y*x)
        """
        xy = self.gp(x, y)
        yx = self.gp(y, x)
        sum_result = xy + yx
        return sum_result * 0.5
    
    def acp(self, x: 'MultiVector', y: 'MultiVector') -> 'MultiVector':
        """Alias for anticommutator_product."""
        return self.anticommutator_product(x, y)
    
    def inversion(self, x: 'MultiVector') -> 'MultiVector':
        """
        Compute the multiplicative inverse of a multivector.
        
        Args:
            x: Multivector to invert
            
        Returns:
            Inverse of x such that x * x.inv() = 1
        """
        # For a multivector x, the inverse is x.reverse() / x.norm_squared()
        # This works for invertible multivectors
        try:
            x_rev = x.reverse()
            norm_sq = self.norm_squared(x)
            
            # Check if norm_squared is zero (non-invertible)
            if isinstance(norm_sq, (int, float, complex)) and abs(norm_sq) < 1e-12:
                raise ZeroDivisionError("Cannot invert multivector with zero norm")
            
            return x_rev / norm_sq
        except (AttributeError, NotImplementedError):
            raise NotImplementedError("Inversion requires reverse() and norm_squared() methods")
    
    def inv(self, x: 'MultiVector') -> 'MultiVector':
        """Alias for inversion."""
        return self.inversion(x)
    
    def division(self, x: 'MultiVector', y: 'MultiVector') -> 'MultiVector':
        """
        Compute the division of two multivectors: x / y = x * y.inv().
        
        Args:
            x: Dividend multivector
            y: Divisor multivector
            
        Returns:
            Division result x / y
        """
        y_inv = self.inversion(y)
        return self.gp(x, y_inv)
    
    def div(self, x: 'MultiVector', y: 'MultiVector') -> 'MultiVector':
        """Alias for division."""
        return self.division(x, y)
    
    def norm_squared(self, x: 'MultiVector') -> Union[float, complex, 'MultiVector']:
        """
        Compute the squared norm of a multivector: x * x.reverse().
        
        Args:
            x: Multivector to compute norm squared for
            
        Returns:
            Squared norm (should be a scalar)
        """
        try:
            x_rev = x.reverse()
            result = self.gp(x, x_rev)
            
            # Extract scalar part if it's a multivector
            if hasattr(result, '_keys') and hasattr(result, '_values'):
                # Find scalar component (key 0)
                if 0 in result._keys:
                    idx = result._keys.index(0)
                    return result._values[idx]
                else:
                    return 0  # No scalar component
            return result
        except (AttributeError, NotImplementedError):
            raise NotImplementedError("norm_squared requires reverse() method")
    
    def normsq(self, x: 'MultiVector') -> Union[float, complex, 'MultiVector']:
        """Alias for norm_squared."""
        return self.norm_squared(x)
    
    def normalization(self, x: 'MultiVector') -> 'MultiVector':
        """
        Normalize a multivector to unit magnitude.
        
        Args:
            x: Multivector to normalize
            
        Returns:
            Normalized multivector x / |x|
        """
        import math
        norm_sq = self.norm_squared(x)
        
        if isinstance(norm_sq, (int, float, complex)):
            if abs(norm_sq) < 1e-12:
                raise ZeroDivisionError("Cannot normalize zero multivector")
            norm = math.sqrt(abs(norm_sq))
            return x / norm
        else:
            # Symbolic case
            import sympy
            norm = sympy.sqrt(norm_sq)
            return x / norm
    
    def sqrt_operation(self, x: 'MultiVector') -> 'MultiVector':
        """
        Compute the square root of a multivector (simplified implementation).
        
        Args:
            x: Multivector to take square root of
            
        Returns:
            Square root of x
        """
        # This is a simplified implementation
        # For a complete implementation, we'd need to handle the general case
        # which can be quite complex for arbitrary multivectors
        
        # For scalars, it's straightforward
        if hasattr(x, '_keys') and len(x._keys) == 1 and x._keys[0] == 0:
            # Pure scalar
            import math
            scalar_val = x._values[0]
            if isinstance(scalar_val, (int, float)) and scalar_val >= 0:
                return self.scalar(math.sqrt(scalar_val))
            else:
                # Complex or symbolic case
                import sympy
                return self.scalar(sympy.sqrt(scalar_val))
        
        # For general multivectors, this is much more complex
        # We'll implement a basic version that works for some cases
        raise NotImplementedError("General multivector square root not yet implemented")
    
    def sqrt(self, x: 'MultiVector') -> 'MultiVector':
        """Alias for sqrt_operation."""
        return self.sqrt_operation(x)
    
    def scalar_product(self, x: 'MultiVector', y: 'MultiVector') -> 'MultiVector':
        """
        Compute the scalar product of two multivectors.
        
        The scalar product extracts only the scalar part of the geometric product.
        It's equivalent to the inner product with diff_func=lambda k: 0.
        
        Args:
            x: First multivector
            y: Second multivector
            
        Returns:
            Scalar product (grade 0 part only)
        """
        from .multivector import MultiVector
        
        if x.algebra != self or y.algebra != self:
            raise ValueError("Multivectors must belong to this algebra")
            
        # Scalar product: only terms where kx ^ ky == 0 (same basis elements)
        result_terms = {}
        
        for kx, vx in zip(x._keys, x._values):
            for ky, vy in zip(y._keys, y._values):
                if kx == ky:  # Same basis element
                    # Get the sign from the Cayley table
                    sign = self.signs[kx, ky] if hasattr(self, 'signs') else 1
                    
                    if sign != 0:
                        key_out = 0  # Always scalar for scalar product
                        
                        # Compute the coefficient
                        if isinstance(vx, (int, float, complex)) and isinstance(vy, (int, float, complex)):
                            coeff = sign * vx * vy
                        else:
                            # Handle symbolic case
                            import sympy
                            coeff = sign * sympy.Mul(vx, vy, evaluate=False) if sign > 0 else -sympy.Mul(vx, vy, evaluate=False)
                        
                        # Accumulate terms
                        if key_out in result_terms:
                            result_terms[key_out] += coeff
                        else:
                            result_terms[key_out] = coeff
        
        # Create result
        if not result_terms:
            return MultiVector.fromkeysvalues(self, keys=(), values=[])
            
        return MultiVector.fromkeysvalues(self, keys=(0,), values=[result_terms[0]])
    
    def sp(self, x: 'MultiVector', y: 'MultiVector') -> 'MultiVector':
        """Alias for scalar_product."""
        return self.scalar_product(x, y)
    
    def projection(self, x: 'MultiVector', y: 'MultiVector') -> 'MultiVector':
        """
        Compute the projection of x onto y: (x · y) * ~y.
        
        Args:
            x: Multivector to project
            y: Multivector to project onto
            
        Returns:
            Projection of x onto y
        """
        # Projection formula: (x | y) * ~y
        inner_prod = self.inner_product(x, y)
        y_reverse = y.reverse()
        return self.gp(inner_prod, y_reverse)
    
    def proj(self, x: 'MultiVector', y: 'MultiVector') -> 'MultiVector':
        """Alias for projection."""
        return self.projection(x, y)
    
    def involute(self, x: 'MultiVector') -> 'MultiVector':
        """
        Compute the grade involution of a multivector.
        
        Grade involution flips the sign of odd grades (1, 3, 5, ...).
        
        Args:
            x: Multivector to apply involution to
            
        Returns:
            Grade involute of x
        """
        from .multivector import MultiVector
        
        if not hasattr(x, '_keys') or not x._keys:
            return x.copy()
        
        # Create new values with involution signs applied
        new_values = []
        for k, v in zip(x._keys, x._values):
            grade = bin(k).count('1')  # Count number of 1s = grade
            
            # Grade involution: flip sign for grades 1, 3, 5, ... (odd grades)
            sign_factor = -1 if grade % 2 == 1 else 1
            
            if isinstance(v, (int, float, complex)):
                new_values.append(sign_factor * v)
            else:
                # Handle symbolic case
                import sympy
                if sign_factor == 1:
                    new_values.append(v)
                else:
                    new_values.append(-v)
        
        return MultiVector.fromkeysvalues(self, x._keys, new_values)
    
    def clifford_conjugate(self, x: 'MultiVector') -> 'MultiVector':
        """
        Compute the Clifford conjugate of a multivector.
        
        Clifford conjugate flips the sign of grades 1 and 2.
        
        Args:
            x: Multivector to apply conjugation to
            
        Returns:
            Clifford conjugate of x
        """
        from .multivector import MultiVector
        
        if not hasattr(x, '_keys') or not x._keys:
            return x.copy()
        
        # Create new values with conjugation signs applied
        new_values = []
        for k, v in zip(x._keys, x._values):
            grade = bin(k).count('1')  # Count number of 1s = grade
            
            # Clifford conjugate: flip sign for grades 1, 2
            sign_factor = -1 if grade in (1, 2) else 1
            
            if isinstance(v, (int, float, complex)):
                new_values.append(sign_factor * v)
            else:
                # Handle symbolic case
                import sympy
                if sign_factor == 1:
                    new_values.append(v)
                else:
                    new_values.append(-v)
        
        return MultiVector.fromkeysvalues(self, x._keys, new_values)
    
    def conjugate(self, x: 'MultiVector') -> 'MultiVector':
        """Alias for clifford_conjugate."""
        return self.clifford_conjugate(x)
    
    def polarity(self, x: 'MultiVector') -> 'MultiVector':
        """
        Compute the polarity (multiplication by pseudoscalar inverse) of a multivector.
        
        Polarity is defined as x * I^(-1) where I is the pseudoscalar.
        
        Args:
            x: Multivector to apply polarity to
            
        Returns:
            Polarity of x
        """
        try:
            pss_inv = self.pss_inv
            return self.gp(x, pss_inv)
        except AttributeError:
            raise NotImplementedError("Polarity requires pseudoscalar inverse (pss_inv)")
    
    def unpolarity(self, x: 'MultiVector') -> 'MultiVector':
        """
        Compute the unpolarity (multiplication by pseudoscalar) of a multivector.
        
        Unpolarity is defined as x * I where I is the pseudoscalar.
        
        Args:
            x: Multivector to apply unpolarity to
            
        Returns:
            Unpolarity of x
        """
        return self.gp(x, self.pss)
    
    def hodge_dual(self, x: 'MultiVector') -> 'MultiVector':
        """
        Compute the Hodge dual of a multivector.
        
        The Hodge dual maps k-vectors to (n-k)-vectors.
        
        Args:
            x: Multivector to compute Hodge dual of
            
        Returns:
            Hodge dual of x
        """
        from .multivector import MultiVector
        
        if not hasattr(x, '_keys') or not x._keys:
            return x.copy()
        
        # Hodge dual: for each basis blade eI, compute eI* = eI * I / |I|^2
        # where I is the pseudoscalar
        result_terms = {}
        algebra_len = len(self)
        key_pss = algebra_len - 1  # All bits set
        
        for k, v in zip(x._keys, x._values):
            # Hodge dual key: complement of k
            key_dual = key_pss ^ k
            
            # Sign from the geometric product with pseudoscalar
            if hasattr(self, 'signs'):
                sign = self.signs.get((k, key_pss), 1)
            else:
                sign = 1
            
            # Apply sign
            if isinstance(v, (int, float, complex)):
                coeff = sign * v
            else:
                # Handle symbolic case
                import sympy
                coeff = sign * v if sign > 0 else -v
            
            if key_dual in result_terms:
                result_terms[key_dual] += coeff
            else:
                result_terms[key_dual] = coeff
        
        # Create result
        if not result_terms:
            return MultiVector.fromkeysvalues(self, keys=(), values=[])
            
        sorted_keys = tuple(sorted(result_terms.keys()))
        sorted_values = [result_terms[k] for k in sorted_keys]
        
        return MultiVector.fromkeysvalues(self, keys=sorted_keys, values=sorted_values)
    
    def hodge(self, x: 'MultiVector') -> 'MultiVector':
        """Alias for hodge_dual."""
        return self.hodge_dual(x)
    
    def unhodge_dual(self, x: 'MultiVector') -> 'MultiVector':
        """
        Compute the inverse Hodge dual of a multivector.
        
        Args:
            x: Multivector to compute inverse Hodge dual of
            
        Returns:
            Inverse Hodge dual of x
        """
        # The inverse Hodge dual is the Hodge dual applied again (in most cases)
        # This is a simplification - the full implementation depends on the metric
        return self.hodge_dual(x)
    
    def unhodge(self, x: 'MultiVector') -> 'MultiVector':
        """Alias for unhodge_dual."""
        return self.unhodge_dual(x)
    
    def outer_exponential(self, x: 'MultiVector') -> 'MultiVector':
        """
        Compute the outer exponential of a multivector: exp_∧(x).
        
        For a k-vector B, exp_∧(B) = 1 + B + B∧B/2! + B∧B∧B/3! + ...
        
        Args:
            x: Multivector to compute outer exponential of
            
        Returns:
            Outer exponential of x
        """
        from .multivector import MultiVector
        import math
        
        # Start with the identity (scalar 1)
        result = self.multivector(values={0: 1})  # Explicit scalar multivector
        
        # For vectors, the series terminates after the first term
        if x.grades == (1,):
            return result + x
        
        # For bivectors and higher, compute the series
        # We'll compute a few terms of the series
        current_term = x.copy()
        factorial = 1
        
        for n in range(1, 8):  # Compute first 7 terms
            factorial *= n
            result = result + current_term / factorial
            
            # Compute next term: current_term ∧ x
            try:
                next_term = self.outer_product(current_term, x)
                # If the outer product is zero, the series terminates
                if next_term.is_zero():
                    break
                current_term = next_term
            except:
                break
        
        return result
    
    def outerexp(self, x: 'MultiVector') -> 'MultiVector':
        """Alias for outer_exponential."""
        return self.outer_exponential(x)
    
    def outer_sine(self, x: 'MultiVector') -> 'MultiVector':
        """
        Compute the outer sine of a multivector.
        
        For a bivector B, sin_∧(B) = B - B∧B∧B/3! + B∧B∧B∧B∧B/5! - ...
        
        Args:
            x: Multivector to compute outer sine of
            
        Returns:
            Outer sine of x
        """
        from .multivector import MultiVector
        import math
        
        # For vectors, sin_∧(v) = v
        if x.grades == (1,):
            return x.copy()
        
        # For bivectors and higher, compute the series
        result = self.scalar(0)
        current_term = x.copy()
        factorial = 1
        sign = 1
        
        for n in range(1, 8, 2):  # Odd terms only: 1, 3, 5, 7, ...
            factorial *= n
            if n > 1:
                factorial *= (n - 1)
            
            result = result + sign * current_term / factorial
            sign *= -1
            
            # Compute next odd term: current_term ∧ x ∧ x
            try:
                next_term = self.outer_product(self.outer_product(current_term, x), x)
                if next_term.is_zero():
                    break
                current_term = next_term
            except:
                break
        
        return result
    
    def outersin(self, x: 'MultiVector') -> 'MultiVector':
        """Alias for outer_sine."""
        return self.outer_sine(x)
    
    def outer_cosine(self, x: 'MultiVector') -> 'MultiVector':
        """
        Compute the outer cosine of a multivector.
        
        For a bivector B, cos_∧(B) = 1 - B∧B/2! + B∧B∧B∧B/4! - ...
        
        Args:
            x: Multivector to compute outer cosine of
            
        Returns:
            Outer cosine of x
        """
        from .multivector import MultiVector
        import math
        
        # For vectors, cos_∧(v) = 1
        if x.grades == (1,):
            return self.multivector(values={0: 1})  # Explicit scalar multivector
        
        # For bivectors and higher, compute the series
        result = self.multivector(values={0: 1})  # Start with scalar 1
        current_term = self.outer_product(x, x)  # Start with x∧x
        factorial = 2
        sign = -1
        
        for n in range(2, 8, 2):  # Even terms: 2, 4, 6, ...
            result = result + sign * current_term / factorial
            sign *= -1
            
            # Compute next even term
            factorial *= (n + 1) * (n + 2)
            
            try:
                next_term = self.outer_product(self.outer_product(current_term, x), x)
                if next_term.is_zero():
                    break
                current_term = next_term
            except:
                break
        
        return result
    
    def outercos(self, x: 'MultiVector') -> 'MultiVector':
        """Alias for outer_cosine."""
        return self.outer_cosine(x)
    
    def outer_tangent(self, x: 'MultiVector') -> 'MultiVector':
        """
        Compute the outer tangent of a multivector: tan_∧(x) = sin_∧(x) / cos_∧(x).
        
        Args:
            x: Multivector to compute outer tangent of
            
        Returns:
            Outer tangent of x
        """
        sin_x = self.outer_sine(x)
        cos_x = self.outer_cosine(x)
        return self.division(sin_x, cos_x)
    
    def outertan(self, x: 'MultiVector') -> 'MultiVector':
        """Alias for outer_tangent."""
        return self.outer_tangent(x)
    
    def add_operation(self, x: 'MultiVector', y: 'MultiVector') -> 'MultiVector':
        """
        Optimized addition operation for multivectors.
        
        This is equivalent to x + y but can be optimized for specific cases.
        
        Args:
            x: First multivector
            y: Second multivector
            
        Returns:
            Sum x + y
        """
        # For now, delegate to the built-in addition
        return x + y
    
    def add(self, x: 'MultiVector', y: 'MultiVector') -> 'MultiVector':
        """Alias for add_operation."""
        return self.add_operation(x, y)
    
    def sub_operation(self, x: 'MultiVector', y: 'MultiVector') -> 'MultiVector':
        """
        Optimized subtraction operation for multivectors.
        
        This is equivalent to x - y but can be optimized for specific cases.
        
        Args:
            x: First multivector
            y: Second multivector
            
        Returns:
            Difference x - y
        """
        # For now, delegate to the built-in subtraction
        return x - y
    
    def sub(self, x: 'MultiVector', y: 'MultiVector') -> 'MultiVector':
        """Alias for sub_operation."""
        return self.sub_operation(x, y)
    
    def neg_operation(self, x: 'MultiVector') -> 'MultiVector':
        """
        Negation operation for multivectors.
        
        Args:
            x: Multivector to negate
            
        Returns:
            Negated multivector -x
        """
        # Delegate to the built-in negation
        return -x
    
    def neg(self, x: 'MultiVector') -> 'MultiVector':
        """Alias for neg_operation."""
        return self.neg_operation(x)

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