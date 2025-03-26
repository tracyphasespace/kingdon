"""
Unified Algebra Module for Geometric Algebra
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
"""

from __future__ import annotations

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

# Field initializer for operator dictionaries
operation_field = partial(field, default_factory=dict, init=False, repr=False, compare=False)


@dataclass
class Algebra:
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
        """Initialize the algebra with computed attributes."""
        # Import here to avoid circular imports at runtime
        from kingdon.multivector import MultiVector
        
        # Set the base namespace for compiled functions
        self.numspace = {'np': np, 'MultiVector': MultiVector}
        
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