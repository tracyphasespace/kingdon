# Kingdon Geometric Algebra Library (Enhanced Fork)

![Version](https://img.shields.io/badge/version-0.2.0-blue)
![Python](https://img.shields.io/badge/python-3.7%2B-blue)
![License](https://img.shields.io/badge/license-MIT-green)

This fork of the Kingdon Geometric Algebra library provides enhanced functionality, improved documentation, and better performance for working with Geometric Algebra in Python.

## Overview

Kingdon is a comprehensive Python implementation of Geometric Algebra (GA), also known as Clifford Algebra. It provides both symbolic and numerical computation capabilities, allowing for coordinate-free manipulation of geometric entities and transformations.

## Key Features

- **Flexible Algebra Creation**: Support for various Clifford algebras including Vector GA, Projective GA, Conformal GA, and custom signatures
- **Multivector Operations**: Complete set of GA operations (geometric product, inner/outer products, rotors, etc.)
- **Performance Optimizations**: Dynamic code generation and common subexpression elimination for efficient numerical computation
- **NumPy Integration**: Enhanced interoperability with NumPy arrays and operations
- **Symbolic Support**: Seamless handling of symbolic expressions alongside numerical values
- **Factory Methods**: Intuitive creation of rotors, translators, reflectors, and other geometric entities

## Improvements in this Fork

This fork includes several key enhancements over the original Kingdon library:

- **Improved Code Organization**: Better structured modules with clear separation of concerns
- **Enhanced Documentation**: More comprehensive docstrings, examples, and type annotations
- **Better NumPy Integration**: Improved handling of array-valued multivectors
- **Computational Robustness**: More robust error handling and edge case management
- **Performance Optimization**: Smarter caching and more efficient implementation of core operations
- **New Convenience Methods**: Added intuitive factory methods for common geometric entities
- **Type Annotations**: Extensive typing for better IDE support and code validation

## Installation

```bash
pip install kingdon-ga
```

## Quick Examples

### Creating an Algebra and Basic Operations

```python
from kingdon.algebra import Algebra

# Create a 3D Geometric Algebra
algebra = Algebra(p=3, q=0, r=0)  # 3D Euclidean space

# Create vectors
v1 = algebra.vector([1, 2, 3])
v2 = algebra.vector([4, 5, 6])

# Geometric product
gp = v1 * v2

# Inner product
inner = v1 | v2

# Outer product
outer = v1 ^ v2

print(f"Geometric product: {gp}")
print(f"Inner product: {inner}")
print(f"Outer product: {outer}")
```

### Creating Rotors and Applying Rotations

```python
import math
from kingdon.algebra import Algebra

# Create a 3D Geometric Algebra
algebra = Algebra(p=3, q=0, r=0)

# Create a vector
v = algebra.vector([1, 0, 0])

# Create a rotor for 90-degree rotation in the xy-plane
R = algebra.rotor(angle=math.pi/2, plane_indices=(1, 2))

# Apply the rotation
v_rotated = R >> v  # Sandwich product: R * v * ~R

print(f"Original vector: {v}")
print(f"Rotor: {R}")
print(f"Rotated vector: {v_rotated}")  # Should be approximately [0, 1, 0]
```

### Projective Geometric Algebra (PGA) Example

```python
from kingdon.algebra import Algebra

# Create a 3D Projective Geometric Algebra (PGA)
pga = Algebra(p=3, q=0, r=1)

# Create a line (bivector in PGA)
line = pga.bivector([1, 0, 0, 0, 0, 1])

# Create a point (weighted vector in PGA)
point = pga.vector([1, 2, 0, 1])

# Compute the meet (intersection)
intersection = line & point

print(f"Line: {line}")
print(f"Point: {point}")
print(f"Intersection: {intersection}")
```

## Applications

The Kingdon library is useful for various applications including:

- **Physics and Engineering**: Relativistic physics, rigid body mechanics, electromagnetic theory
- **Computer Graphics**: Camera transformations, projections, and rigid body dynamics
- **Robotics**: Kinematics, path planning, and spatial reasoning
- **Machine Learning**: Geometric feature representation and manifold learning

## Contributing

Contributions to the library are welcome! Please feel free to submit issues and pull requests.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

This fork builds upon the excellent work of the original Kingdon library creators, enhancing and extending their foundation to provide a more robust and user-friendly implementation of Geometric Algebra in Python.
