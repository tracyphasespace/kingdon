"""
Factory Functions for Geometric Algebra Objects
=================================================

This module provides a simple, high-level interface for creating common
Geometric Algebra objects such as vectors, bivectors, and scalars. These
factory functions help reduce the complexity of directly calling the
underlying multivector constructors.

The module includes functions for creating:
- Vectors (grade 1)
- Bivectors (grade 2)
- Trivectors (grade 3)
- Scalars (grade 0)
- General multivectors with specific grades
- Special elements like pseudoscalars, rotors, and reflectors

Example:
    >>> from kingdon.algebra import Algebra
    >>> from kingdon.ga_factory import create_vector, create_bivector
    >>> alg = Algebra(p=3, q=0, r=1)  # 3D PGA
    >>>
    >>> # Create a vector
    >>> v = create_vector(alg, values=[1, 2, 3, 0], name="v")
    >>>
    >>> # Create a bivector
    >>> B = create_bivector(alg, values=[0, 1, 0, 0, 0, 0], name="B")
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union
import numpy as np



def create_vector(algebra: Any, 
                  values: Union[Sequence[Any], np.ndarray, Dict[int, Any]], 
                  name: Optional[str] = None) -> Any:
    """
    Create a vector (grade 1 multivector) in the given algebra.

    Args:
        algebra: The Geometric Algebra instance
        values: The coefficients for the vector's components
        name: Optional name for symbolic representation

    Returns:
        A multivector representing a vector

    Example:
        >>> from kingdon.algebra import Algebra
        >>> from kingdon.ga_factory import create_vector
        >>> alg = Algebra(p=3, q=0, r=1)  # 3D PGA
        >>> v = create_vector(alg, values=[1, 0, 0, 0], name="v")
        >>> print(v)
        v = 1𝐞₀
    """
    # Use algebra's vector method directly to avoid circular imports
    return algebra.vector(values, name=name)


def create_bivector(algebra: Any, 
                   values: Union[Sequence[Any], np.ndarray, Dict[int, Any]], 
                   name: Optional[str] = None) -> Any:
    """
    Create a bivector (grade 2 multivector) in the given algebra.

    Args:
        algebra: The Geometric Algebra instance
        values: The coefficients for the bivector's components
        name: Optional name for symbolic representation

    Returns:
        A multivector representing a bivector

    Example:
        >>> from kingdon.algebra import Algebra
        >>> from kingdon.ga_factory import create_bivector
        >>> alg = Algebra(p=3, q=0, r=0)  # 3D Euclidean GA
        >>> # Create bivector representing rotation in xy-plane
        >>> B = create_bivector(alg, values=[1, 0, 0], name="B")
        >>> print(B)
        B = 1𝐞₁₂
    """
    # Use algebra's multivector method directly with grade 2
    return algebra.multivector(values=values, grades=(2,), name=name)


def create_trivector(algebra: Any, 
                    values: Union[Sequence[Any], np.ndarray, Dict[int, Any]], 
                    name: Optional[str] = None) -> Any:
    """
    Create a trivector (grade 3 multivector) in the given algebra.

    Args:
        algebra: The Geometric Algebra instance
        values: The coefficients for the trivector's components
        name: Optional name for symbolic representation

    Returns:
        A multivector representing a trivector

    Example:
        >>> from kingdon.algebra import Algebra
        >>> from kingdon.ga_factory import create_trivector
        >>> alg = Algebra(p=3, q=0, r=1)  # 3D PGA
        >>> # Create a trivector (pseudovector in 3D)
        >>> I3 = create_trivector(alg, values=[1], name="I3")
        >>> print(I3)
        I3 = 1𝐞₁₂₃
    """
    # Use algebra's multivector method directly with grade 3
    return algebra.multivector(values=values, grades=(3,), name=name)


def create_scalar(algebra: Any, value: Any, name: Optional[str] = None) -> Any:
    """
    Create a scalar (grade 0 multivector) in the given algebra.

    Args:
        algebra: The Geometric Algebra instance
        value: The scalar value (single coefficient)
        name: Optional name for symbolic representation

    Returns:
        A multivector representing a scalar

    Example:
        >>> from kingdon.algebra import Algebra
        >>> from kingdon.ga_factory import create_scalar
        >>> alg = Algebra(p=3, q=0, r=1)
        >>> s = create_scalar(alg, value=5, name="s")
        >>> print(s)
        s = 5
    """
    # Use algebra's scalar method directly
    return algebra.scalar(value, name=name)


def create_multivector(algebra: Any, 
                      grades: Tuple[int, ...], 
                      values: Union[Sequence[Any], np.ndarray, Dict[str, Any]], 
                      name: Optional[str] = None) -> Any:
    """
    Create a multivector with specific grades in the given algebra.

    Args:
        algebra: The Geometric Algebra instance
        grades: Tuple of grades to include in the multivector
        values: The coefficients for the multivector components
        name: Optional name for symbolic representation

    Returns:
        A multivector with the specified grades

    Example:
        >>> from kingdon.algebra import Algebra
        >>> from kingdon.ga_factory import create_multivector
        >>> alg = Algebra(p=3, q=0, r=0)  # 3D Euclidean GA
        >>> # Create a multivector with scalar and vector parts
        >>> mv = create_multivector(alg, grades=(0, 1), values=[1, 2, 3, 4], name="mv")
        >>> print(mv)
        mv = 1 + 2𝐞₁ + 3𝐞₂ + 4𝐞₃
    """
    # Use algebra's multivector method directly
    return algebra.multivector(values=values, grades=grades, name=name)


def create_blade(algebra: Any, 
                indices: Sequence[int], 
                value: Any = 1, 
                name: Optional[str] = None) -> Any:
    """
    Create a specific basis blade by indices.

    Args:
        algebra: The Geometric Algebra instance
        indices: The indices of the basis vectors to wedge together
            - For example, [1, 2] creates the e₁₂ blade
        value: Coefficient for the blade (default is 1)
        name: Optional name for symbolic representation

    Returns:
        A multivector representing the specified basis blade

    Example:
        >>> from kingdon.algebra import Algebra
        >>> from kingdon.ga_factory import create_blade
        >>> alg = Algebra(p=3, q=0, r=0)  # 3D Euclidean GA
        >>> # Create the e₁₃ blade with coefficient 2
        >>> B = create_blade(alg, indices=[1, 3], value=2, name="B")
        >>> print(B)
        B = 2𝐞₁₃
    """
    if not indices:
        return create_scalar(algebra, value, name)
    
    # Try to use the get_blade_by_indices method if available
    if hasattr(algebra, 'get_blade_by_indices'):
        try:
            blade = algebra.get_blade_by_indices(*indices)
            return blade * value
        except (AttributeError, TypeError):
            pass
    
    # Try using the blades attribute if available
    if hasattr(algebra, 'blades'):
        # Create blade name (e.g., "e12" for e₁₂)
        sorted_indices = sorted(indices)
        blade_name = f"e{''.join(str(i) for i in sorted_indices)}"
        
        try:
            # Try attribute access first
            blade = getattr(algebra.blades, blade_name)
            return blade * value
        except AttributeError:
            try:
                # Try dictionary access
                blade = algebra.blades[blade_name]
                return blade * value
            except (KeyError, TypeError):
                pass
    
    # As a fallback, create blade using vector operations
    if len(indices) == 1:
        # For grade 1, just create a vector with the appropriate component
        idx = indices[0]
        vec_values = [0] * algebra.d
        adjusted_idx = idx - algebra.start_index
        if 0 <= adjusted_idx < algebra.d:
            vec_values[adjusted_idx] = value
        return create_vector(algebra, values=vec_values, name=name)
    else:
        # For higher grades, compose using outer product
        e_vectors = []
        for idx in indices:
            vec_values = [0] * algebra.d
            adjusted_idx = idx - algebra.start_index
            if 0 <= adjusted_idx < algebra.d:
                vec_values[adjusted_idx] = 1
            e_vectors.append(create_vector(algebra, values=vec_values))
        
        # Use the outer product to combine the vectors
        result = e_vectors[0]
        for vec in e_vectors[1:]:
            result = result ^ vec
        
        return result * value


def create_pseudoscalar(algebra: Any, 
                       value: Any = 1, 
                       name: Optional[str] = None) -> Any:
    """
    Create the pseudoscalar (highest grade element) for the given algebra.

    The pseudoscalar represents the oriented volume element of the space.

    Args:
        algebra: The Geometric Algebra instance
        value: Coefficient for the pseudoscalar (default is 1)
        name: Optional name for symbolic representation

    Returns:
        A multivector representing the pseudoscalar

    Example:
        >>> from kingdon.algebra import Algebra
        >>> from kingdon.ga_factory import create_pseudoscalar
        >>> alg = Algebra(p=3, q=0, r=0)  # 3D Euclidean GA
        >>> I = create_pseudoscalar(alg, name="I")
        >>> print(I)
        I = 1𝐞₁₂₃
    """
    # Try to access the pseudoscalar directly if available
    if hasattr(algebra, 'pss'):
        pss = algebra.pss
        if value != 1:
            pss = pss * value
        
        # Set name if provided
        if name and hasattr(pss, 'name'):
            pss.name = name
        
        return pss
    
    # Otherwise create using the highest grade
    d = algebra.d
    return algebra.multivector(values=[value], grades=(d,), name=name)


def create_rotor(algebra: Any, 
                angle: float, 
                plane_indices: Tuple[int, int], 
                name: Optional[str] = None) -> Any:
    """
    Create a rotor representing a rotation in a specific plane.

    A rotor is of the form R = cos(θ/2) + sin(θ/2)Bᵢⱼ where Bᵢⱼ is a unit bivector.

    Args:
        algebra: The Geometric Algebra instance
        angle: Rotation angle in radians
        plane_indices: Tuple of two indices defining the rotation plane
            - For example, (1, 2) for rotation in the e₁e₂ plane
        name: Optional name for symbolic representation

    Returns:
        A multivector representing the rotor

    Example:
        >>> import math
        >>> from kingdon.algebra import Algebra
        >>> from kingdon.ga_factory import create_rotor
        >>> alg = Algebra(p=3, q=0, r=0)  # 3D Euclidean GA
        >>> # Create a rotor for 90° rotation in the xy-plane
        >>> R = create_rotor(alg, angle=math.pi/2, plane_indices=(1, 2), name="R")
        >>> print(R)
        R = 0.7071067811865476 + 0.7071067811865475𝐞₁₂
    """
    # Create scalar part cos(θ/2)
    scalar = create_scalar(algebra, value=math.cos(angle/2))
    
    # Create bivector part sin(θ/2)Bᵢⱼ
    i, j = sorted(plane_indices)
    
    # Generate the bivector key 
    # e.g. for (1,2) we need e12 which corresponds to binary key 3
    binary_key = (1 << (i-1)) | (1 << (j-1))
    bivector = create_blade(algebra, indices=[i, j], value=math.sin(angle/2))
    
    # Combine scalar and bivector parts
    from kingdon.multivector import MultiVector
    if isinstance(scalar, MultiVector) and isinstance(bivector, MultiVector):
        # Get the scalar and bivector values
        scalar_value = scalar._values[0]
        bivector_key = bivector._keys[0]
        bivector_value = bivector._values[0]
        
        # Create a new multivector with both components
        rotor = MultiVector.fromkeysvalues(
            algebra,
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

def create_translator(algebra: Any, 
                     direction: Sequence[float], 
                     distance: float = 1.0, 
                     name: Optional[str] = None) -> Any:
    """
    Create a translator in Projective Geometric Algebra (PGA).

    A translator is of the form T = 1 + (d/2)e₀ᵢvᵢ where e₀ᵢ is a null vector direction.
    This only works in algebras with at least one null dimension (r ≥ 1).

    Args:
        algebra: The Geometric Algebra instance (must have r ≥ 1)
        direction: Direction vector components (excluding null component)
        distance: Translation distance
        name: Optional name for symbolic representation

    Returns:
        A multivector representing the translator

    Example:
        >>> from kingdon.algebra import Algebra
        >>> from kingdon.ga_factory import create_translator
        >>> alg = Algebra(p=3, q=0, r=1)  # 3D PGA
        >>> # Create a translator for 2 units in the x direction
        >>> T = create_translator(alg, direction=[1, 0, 0], distance=2, name="T")
        >>> print(T)
        T = 1 + 1.0𝐞₀₁
    """
    # Validate that the algebra has at least one null dimension
    if hasattr(algebra, 'r') and algebra.r < 1:
        raise ValueError("Translator requires an algebra with at least one null dimension (r ≥ 1)")
    
    # Create scalar part (1)
    translator = create_scalar(algebra, value=1)
    
    # Create bivector parts (d/2)e₀ᵢvᵢ
    for i, component in enumerate(direction, start=1):
        if component == 0:
            continue
            
        # Create bivector e₀ᵢ with coefficient (d/2)*vᵢ
        blade = create_blade(algebra, indices=[0, i], value=(distance/2) * component)
        translator = translator + blade
    
    # Set name if needed
    if name and hasattr(translator, 'name'):
        translator.name = name
        
    return translator


import math
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union
import numpy as np

def create_reflector(algebra: Any, 
                    normal: Sequence[float], 
                    name: Optional[str] = None) -> Any:
    """
    Create a reflector (versor representing reflection in a hyperplane).

    A reflector is simply a normalized vector perpendicular to the reflection hyperplane.

    Args:
        algebra: The Geometric Algebra instance
        normal: Normal vector components defining the reflection hyperplane
        name: Optional name for symbolic representation

    Returns:
        A multivector representing the reflector

    Example:
        >>> from kingdon.algebra import Algebra
        >>> from kingdon.ga_factory import create_reflector
        >>> alg = Algebra(p=3, q=0, r=0)  # 3D Euclidean GA
        >>> # Create a reflector for the xy-plane (normal = z-axis)
        >>> R = create_reflector(alg, normal=[0, 0, 1], name="R")
        >>> print(R)
        R = 1.0𝐞₃
    """
    # Validate input - but be more flexible with dimension
    if len(normal) < algebra.d - 1:
        raise ValueError(f"Normal vector must have at least {algebra.d-1} components, got {len(normal)}")
    
    # For PGA (r=1), we need to handle the special case where the normal might not include the homogeneous coordinate
    if hasattr(algebra, 'r') and algebra.r == 1 and len(normal) == algebra.d - 1:
        # Add a 0 for the homogeneous coordinate
        normal = list(normal) + [0]
    
    # Pad with zeros if needed
    if len(normal) < algebra.d:
        normal = list(normal) + [0] * (algebra.d - len(normal))
    
    # Check if the normal vector is zero
    norm_sq = sum(x*x for x in normal)
    if norm_sq == 0:
        raise ValueError("Normal vector cannot be zero")
    
    # Create the normal vector
    vector = create_vector(algebra, values=normal, name=name)
    
    # Calculate the norm manually for normalization
    norm = math.sqrt(norm_sq)
    
    # Normalize by manually dividing - this avoids needing the norm() method
    from kingdon.multivector import MultiVector
    if isinstance(vector, MultiVector):
        # Create a new multivector with normalized values
        normalized_values = [v / norm for v in vector._values]
        return MultiVector.fromkeysvalues(algebra, vector._keys, normalized_values)
    else:
        # Fallback to division
        return vector / norm


def create_motor(algebra: Any, 
                rotation: Optional[Tuple[float, Tuple[int, int]]] = None,
                translation: Optional[Tuple[Sequence[float], float]] = None,
                name: Optional[str] = None) -> Any:
    """
    Create a motor (combined rotation and translation) in PGA.

    A motor is the product of a rotor and a translator and represents a rigid body motion.

    Args:
        algebra: The Geometric Algebra instance (must have r ≥ 1 for translation)
        rotation: Optional tuple (angle, plane_indices) for rotation component
        translation: Optional tuple (direction, distance) for translation component
        name: Optional name for symbolic representation

    Returns:
        A multivector representing the motor

    Example:
        >>> import math
        >>> from kingdon.algebra import Algebra
        >>> from kingdon.ga_factory import create_motor
        >>> alg = Algebra(p=3, q=0, r=1)  # 3D PGA
        >>> # Create a motor for 90° rotation in xy-plane followed by translation in x
        >>> M = create_motor(
        ...     alg, 
        ...     rotation=(math.pi/2, (1, 2)),
        ...     translation=([1, 0, 0], 2),
        ...     name="M"
        ... )
    """
    # Check if either rotation or translation is provided
    if not rotation and not translation:
        return create_scalar(algebra, value=1, name=name)
    
    # Check if the algebra has a null dimension for translation
    if translation and hasattr(algebra, 'r') and algebra.r < 1:
        raise ValueError("Translation requires an algebra with at least one null dimension (r ≥ 1)")
    
    # Start with identity motor
    motor = create_scalar(algebra, value=1)
    
    # Apply rotation if specified
    if rotation:
        angle, plane_indices = rotation
        rotor = create_rotor(algebra, angle, plane_indices)
        motor = motor * rotor
    
    # Apply translation if specified
    if translation:
        direction, distance = translation
        translator = create_translator(algebra, direction, distance)
        motor = motor * translator
    
    # Set name if needed
    if name and hasattr(motor, 'name'):
        motor.name = name
        
    return motor