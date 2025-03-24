"""
Base Algebra Module for Geometric Algebra Core Functionality

This module provides the foundational Algebra class with core initialization 
and fundamental geometric algebra operations.
"""

from __future__ import annotations

import operator
import re
from itertools import product
from functools import partial, reduce, lru_cache
from collections import Counter
from dataclasses import dataclass, field

from typing import List, Dict, Optional, Any, Tuple, Callable, Union, Set, Type, TypeVar, cast

import numpy as np
import sympy

# Import code generation functions and supporting classes
from kingdon.codegen import (
    codegen_gp, codegen_sw, codegen_cp, codegen_ip, codegen_op, codegen_div,
    codegen_rp, codegen_acp, codegen_proj, codegen_sp, codegen_lc, codegen_inv,
    codegen_rc, codegen_normsq, codegen_add, codegen_sub, codegen_neg, codegen_reverse,
    codegen_involute, codegen_conjugate, codegen_sqrt, codegen_outerexp, codegen_outersin,
    codegen_outercos, codegen_outertan, codegen_polarity, codegen_unpolarity,
    codegen_hodge, codegen_unhodge
)
from kingdon.operator_dict import OperatorDict, UnaryOperatorDict, BladeDict

# Type variable for operator dictionaries
T = TypeVar('T', bound=OperatorDict)

# Field initializer for operator dictionaries
operation_field = partial(field, default_factory=dict, init=False, repr=False, compare=False)

@dataclass
class AlgebraBase:
    """
    Base class for Geometric (Clifford) Algebra core functionality.

    This class provides the fundamental structure and initialization 
    for geometric algebras, handling signature, dimensions, and basic operations.
    """
    p: int = 0
    q: int = 0
    r: int = 0
    signature: Optional[List[int]] = None
    start_index: Optional[int] = None
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
        # Signature and dimension setup methods
        self._setup_signature()
        self._setup_basis_mappings()
        self._precompute_grade_indices()
        self._compute_multiplication_signs()
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
        """Initialize the blades dictionary and create the pseudoscalar."""
        # Create a dictionary of all basis blades
        self.blades = BladeDict(self, lazy=True)
        
        # Create pseudoscalar (highest grade element)
        pseudoscalar_name = f"e{''.join(str(i + self.start_index) for i in range(self.d))}"
        try:
            self.pss = self.blades[pseudoscalar_name]
        except Exception as e:
            print(f"Warning: Could not create pseudoscalar {pseudoscalar_name}: {e}")
            # Create a fallback pseudoscalar with the highest binary key
            highest_key = (1 << self.d) - 1
            if hasattr(self, 'multivector'):
                self.pss = self.multivector(values=[1], keys=(highest_key,))
            else:
                from kingdon.multivector import MultiVector
                mv = object.__new__(MultiVector)
                mv.algebra = self
                mv._keys = (highest_key,)
                mv._values = [1]
                mv.shape = ()
                self.pss = mv

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