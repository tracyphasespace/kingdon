# Suggested filename: test_functional_api.py
# Refactored from Testing_Kingdon_changes.py

import pytest
import math
import numpy as np
import sys
import os

# Assuming kingdon is installed or available in the python path
from kingdon import Algebra, MultiVector, __version__
from kingdon.algebra import PGA3D, CGA3D, VGA3D # Convenience classes
from kingdon.ga_factory import (
    create_vector,
    create_bivector,
    create_scalar,
    create_blade,
    create_pseudoscalar,
    create_rotor,
    create_translator,
    create_reflector,
    create_motor
)

# --- Fixtures ---

@pytest.fixture(scope='module')
def vga3d():
    """Provides a 3D Vector Geometric Algebra R(3,0,0) instance."""
    return Algebra(p=3, q=0, r=0)

@pytest.fixture(scope='module')
def pga3d():
    """Provides a 3D Projective Geometric Algebra R(3,0,1) instance."""
    return Algebra(p=3, q=0, r=1)

@pytest.fixture(scope='module')
def cga3d():
    """Provides a 3D Conformal Geometric Algebra R(4,1,0) instance."""
    return Algebra(p=4, q=1, r=0)

@pytest.fixture(scope='module')
def sta():
    """Provides a Spacetime Algebra R(1,3,0) instance."""
    return Algebra(p=1, q=3, r=0)

# --- Tests ---

def test_library_import_and_version():
    """Tests basic import of the kingdon library and checks version."""
    import kingdon
    assert __version__ is not None, "Kingdon version should be accessible."
    print(f"Kingdon version: {__version__}") # Keep informative print

def test_module_imports():
    """Tests imports of key modules and classes."""
    from kingdon.algebra import Algebra
    from kingdon.multivector import MultiVector
    import kingdon.ga_factory
    # If imports succeed, the test passes implicitly. Add assertions if needed.
    assert True

def test_create_algebras(vga3d, pga3d, cga3d, sta):
    """Tests creation of various standard Geometric Algebras via API."""
    # Check properties of fixture algebras
    assert vga3d.signature.tolist() == [1, 1, 1]
    assert pga3d.signature.tolist() == [1, 1, 1, 0]
    assert cga3d.signature.tolist() == [1, 1, 1, 1, -1]
    assert sta.signature.tolist() == [1, -1, -1, -1]

    # Test convenience functions (already imported)
    pga_conv = PGA3D()
    assert pga_conv.signature.tolist() == pga3d.signature.tolist()

    cga_conv = CGA3D()
    assert cga_conv.signature.tolist() == cga3d.signature.tolist()

    vga_conv = VGA3D()
    assert vga_conv.signature.tolist() == vga3d.signature.tolist()

    print(f"Algebras created: VGA3D={vga3d}, PGA3D={pga3d}, CGA3D={cga3d}, STA={sta}") # Informative

def test_create_multivectors(vga3d):
    """Tests creation of various multivector types via API."""
    alg = vga3d

    # Create vectors using blades attribute API
    v1 = alg.blades.e1
    v2 = alg.blades.e2
    v3 = alg.blades.e3
    assert v1.grades == (1,) and v1.e1 == 1
    assert v2.grades == (1,) and v2.e2 == 1
    assert v3.grades == (1,) and v3.e3 == 1

    # Create bivector using generic multivector method API
    # Specify components for e12, e13, e23
    biv = alg.multivector({3: 1, 5: 2, 6: 3}) # keys 3=e12, 5=e13, 6=e23
    assert biv.grades == (2,)
    assert biv.e12 == 1 and biv.e13 == 2 and biv.e23 == 3

    # Create bivector using specific bivector method API
    biv2 = alg.bivector([1, 2, 3]) # values for e12, e13, e23
    assert biv == biv2 # Should be equal

    # Create scalar using API
    s = alg.scalar(5)
    assert s.grades == (0,) and s.e == 5
    assert s == 5 # Test equality with number

    # Create blade using factory function API
    # Note: Factory functions might be deprecated if alg.blade is preferred
    e23_factory = create_blade(alg, indices=(2, 3), value=2)
    e23_method = alg.blade((2, 3), value=2)
    e23_direct = alg.blades.e23 * 2
    assert e23_factory.grades == (2,) and e23_factory.e23 == 2
    assert e23_method == e23_factory
    assert e23_direct == e23_factory

    # Access blades from the blade dictionary API
    e13 = alg.blades.e13
    assert e13.grades == (2,) and e13.e13 == 1

def test_basic_operations(vga3d):
    """Tests basic multivector operations using the API."""
    alg = vga3d
    v1, v2, v3 = alg.blades.e1, alg.blades.e2, alg.blades.e3
    s = alg.scalar(5)

    # Geometric product (*) API
    assert v1 * v1 == 1, "v1 * v1 should be scalar 1 in R(3)"
    assert v1 * v2 == alg.blades.e12, "v1 * v2 should be e12 in R(3)"
    assert v1 * v2 * v3 == alg.blades.e123, "v1 * v2 * v3 should be e123"

    # Outer product (^) API
    assert v1 ^ v2 == alg.blades.e12, "v1 ^ v2 should be e12"
    assert v1 ^ v1 == 0, "v1 ^ v1 should be 0"
    assert v1 ^ v2 ^ v3 == alg.blades.e123, "v1 ^ v2 ^ v3 should be e123"
    assert v1 ^ (v1^v2) == 0, "v1 ^ (v1^v2) should be 0"

    # Inner product (|) API (Left contraction assumed)
    assert v1 | v1 == 1, "v1 | v1 should be 1"
    assert v1 | v2 == 0, "v1 | v2 should be 0"
    assert v1 | (v1^v2) == v2, "v1 | (v1^v2) should be v2" # v1 | e12 = e2
    assert (v1^v2) | v1 == -v2, "(v1^v2) | v1 should be -v2" # e12 | e1 = -e2 (Check convention)

    # Addition (+) API
    sum_v = v1 + v2
    assert sum_v.grades == (1,) and sum_v.e1 == 1 and sum_v.e2 == 1

    # Subtraction (-) API
    diff_v = v1 - v2
    assert diff_v.grades == (1,) and diff_v.e1 == 1 and diff_v.e2 == -1

    # Scalar multiplication (*) API
    scaled_v = s * v1
    assert scaled_v.grades == (1,) and scaled_v.e1 == 5
    scaled_v2 = v1 * s
    assert scaled_v2 == scaled_v # Commutativity with scalar

    # Reverse (~) API
    assert ~v1 == v1, "~v1 should be v1"
    assert ~alg.blades.e12 == -alg.blades.e12, "~e12 should be -e12"
    assert ~alg.blades.e123 == -alg.blades.e123, "~e123 should be -e123" # Grade 3: k(k-1)/2 = 3*2/2 = 3 (odd) -> sign flip? Check definition. Reverse sign depends on grade mod 4. Gr 3: (-1)**(3*2/2) = -1. Correct.
    assert ~(s + v1 + alg.blades.e12) == s + v1 - alg.blades.e12

    # Test division (/) API
    assert v1 / v1 == 1
    assert v1 / s == v1 * (1/5) == alg.scalar(0.2) * v1
    assert alg.blades.e12 / v2 == alg.blades.e12 * (~v2 / (v2*v2)) == alg.blades.e12 * v2 == -v1 # e12 * e2 = -e1


def test_grades_and_projections(vga3d):
    """Tests grade extraction and projection operation API."""
    alg = vga3d
    mv = 1 + 2*alg.blades.e1 + 3*alg.blades.e12 + 4*alg.blades.e123

    # Test grades property API
    assert set(mv.grades) == {0, 1, 2, 3}

    # Test grade() method API
    assert mv.grade(0) == 1
    assert mv.grade(1) == 2*alg.blades.e1
    assert mv.grade(2) == 3*alg.blades.e12
    assert mv.grade(3) == 4*alg.blades.e123
    assert mv.grade(4) == 0 # Check non-existent grade

    # Test convenience methods API
    assert mv.scalar() == 1
    assert mv.vector() == 2*alg.blades.e1
    assert mv.bivector() == 3*alg.blades.e12
    # Add trivector or pseudoscalar method test if available
    # assert mv.trivector() == 4*alg.blades.e123

def test_specialized_factories(pga3d):
    """Tests specialized factory functions API for GA objects."""
    # Use PGA for testing these factories
    pga = pga3d

    # Test pseudoscalar creation API
    ps = pga.pseudoscalar() # Assuming algebra method exists
    # Or use factory: ps_factory = create_pseudoscalar(pga)
    # PGA3D pseudoscalar is grade 4 (e0123) if e0 is last basis vector
    # Check dimension and grades
    assert ps.grades == (pga.d,), f"Pseudoscalar should be grade {pga.d}, got {ps.grades}"

    # Test rotor creation API
    # Rotor is exp(angle/2 * plane), grades (0, 2)
    rotor = create_rotor(pga, angle=math.pi/8, plane_indices=(1, 2))
    assert rotor.grades == (0, 2), f"Rotor should have grades 0 and 2, got {rotor.grades}"
    # Check normalization R*~R = 1
    assert abs((rotor * ~rotor).e - 1) < 1e-12, "Rotor should be normalized"

    # Test translator creation API
    # Translator T = 1 + d/2 * ideal_line, grades (0, 2)
    translator = create_translator(pga, direction=[1, 0, 0]) # Direction maps to ideal line e10?
    assert translator.grades == (0, 2), f"Translator should have grades 0 and 2, got {translator.grades}"
    # Check normalization T*~T = 1
    assert abs((translator * ~translator).e - 1) < 1e-12, "Translator should be normalized"

    # Test reflector creation API
    # Reflector is a plane (vector in dual?) normal vector, grade 1
    reflector = create_reflector(pga, normal=[0, 0, 1, 0]) # Plane x3=0 ? Should be e3?
    # Factory likely creates the plane directly. Assuming normal=[a,b,c,d] -> a*e1+b*e2+c*e3+d*e0
    # For PGA, reflector is often represented by the plane itself.
    plane = pga.vector([0, 0, 1, 0]) # Create plane e3 via API
    reflector_from_plane = create_reflector(pga, plane=plane)
    assert reflector_from_plane.grades == (1,), f"Reflector plane should be grade 1, got {reflector.grades}"
    assert reflector_from_plane == plane
    # Check normalization P*P = 1 (for Euclidean planes)
    assert plane*plane == 1, "Normalized plane reflection P*P should be 1" # Only if P is normalized e3

    # Test motor creation API (Rotor * Translator)
    motor = create_motor(pga,
                        rotation=(math.pi/20, (1, 2)), # Angle, plane indices
                        translation=([1, 0, 0], 0.5)) # Direction, distance
    # Motor has grades (0, 2)
    assert motor.grades == (0, 2), f"Motor should have grades 0 and 2, got {motor.grades}"
    # Check normalization M*~M = 1
    assert abs((motor * ~motor).e - 1) < 1e-12, "Motor should be normalized"


def test_pga_operations(pga3d):
    """Tests operations specific to Projective Geometric Algebra API."""
    pga = pga3d
    e1, e2, e3, e0 = pga.blades.e1, pga.blades.e2, pga.blades.e3, pga.blades.e0

    # Create a point (vector normalized with e0=1) using API
    point = e1 + 2*e2 + 3*e3 + e0 # pga.vector([1, 2, 3, 1])
    assert point.grades == (1,)

    # Create a plane (vector with e0=0 for ideal, or non-zero for real?)
    # Plane ax+by+cz+d=0 -> dual is a*e1+b*e2+c*e3+d*e0
    plane = e3 + e0 # Plane x3 + 1 = 0, use pga.vector([0,0,1,1])
    assert plane.grades == (1,)

    # Meet operation (intersection) using Regressive Product (&) API
    # point & plane = should be a line (grade 2)
    intersection_line = point & plane
    assert intersection_line.grades == (2,)

    # Join operation (span) using Outer Product (^) API
    # point ^ plane = should be the volume/pseudoscalar (grade 3 in dual? grade 3 here?)
    # PGA3D: point=vec(gr1), plane=vec(gr1). Join is point ^ plane (gr2)
    # Let's join two points
    point2 = -e1 + e2 - e3 + e0
    line_join = point ^ point2
    assert line_join.grades == (2,)

def test_array_valued_functional(vga3d):
    """Tests functional aspects of array-valued multivectors using API."""
    alg = vga3d
    e1 = alg.blades.e1

    # Create array-valued MV using API (scalar * basis_vector)
    identity = np.eye(3) # 3x3 identity matrix
    # We want a MV where each component is e1, but coeffs form identity matrix
    # Need to create this carefully using the API.
    # Example: 3x3 array of scalars times e1
    mv_array = e1 * identity
    assert isinstance(mv_array, MultiVector)
    assert mv_array.shape == (3, 3)

    # Test application of NumPy ufuncs (API)
    mv_squared = np.square(mv_array) # Square the coefficients element-wise
    assert isinstance(mv_squared, MultiVector)
    assert mv_squared.shape == (3, 3)
    # Check a value: element [0,0] should be (1*e1)^2 -> API should handle this?
    # np.square acts on coeffs: (1*e1) -> 1^2 * e1 = 1*e1 ?
    # Ufunc likely applies to coefficients only.
    expected_coeffs_squared = np.square(identity)
    # Check value of element [0,0] after ufunc
    assert mv_squared[0, 0] == e1 * expected_coeffs_squared[0, 0] # 1*e1
    # Check value of element [1,1] after ufunc
    assert mv_squared[1, 1] == e1 * expected_coeffs_squared[1, 1] # 1*e1
    # Check value of element [0,1] after ufunc
    assert mv_squared[0, 1] == e1 * expected_coeffs_squared[0, 1] # 0*e1 = 0


def test_error_recovery(vga3d):
    """Tests error handling for invalid API usage."""
    alg = vga3d

    # Intentionally create an error (different length keys and values) via API
    with pytest.raises((ValueError, TypeError)): # Expect ValueError or TypeError
        # Use generic multivector constructor API
        invalid_mv = alg.multivector(keys=(1, 2), values=[1])

    # Verify we can still create valid objects after an error
    valid_mv = alg.multivector(keys=(1,), values=[1])
    assert valid_mv == alg.blades.e1 # Check result via API