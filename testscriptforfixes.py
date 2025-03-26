"""
Test script to verify the fixes for the Kingdon Geometric Algebra library.
"""

import unittest
import numpy as np
from kingdon.algebra import Algebra
from kingdon.multivector import MultiVector

class TestKingdonFixes(unittest.TestCase):
    def setUp(self):
        """Set up test environment."""
        self.vga3d = Algebra(p=3, q=0, r=0)  # 3D Euclidean
        self.pga3d = Algebra(p=3, q=0, r=1)  # 3D Projective
        
        # Create a factory function for multivectors
        self.mk = lambda **kwargs: MultiVector.fromkeysvalues(
            self.vga3d, 
            keys=tuple(k for k in kwargs.keys() if k != "e"),
            values=[v for k, v in kwargs.items() if k != "e"]
        ) if "e" not in kwargs else MultiVector.fromkeysvalues(
            self.vga3d, keys=(0,), values=[kwargs["e"]]
        )
    
    def test_conjugate(self):
        """Test the Clifford conjugate operation."""
        # Create a scalar
        scalar = self.mk(e=2)
        # Conjugate should not change scalars
        conj_scalar = scalar.conjugate()
        self.assertEqual(conj_scalar.e, 2)
        
        # Create a vector 
        vector = self.mk(**{1: 1, 2: 2, 4: 3})
        # Conjugate should negate vectors
        conj_vector = vector.conjugate()
        self.assertEqual(conj_vector.get(1), -1)
        self.assertEqual(conj_vector.get(2), -2)
        self.assertEqual(conj_vector.get(4), -3)
        
        # Create a bivector
        bivector = self.mk(**{3: 1, 5: 2, 6: 3})
        # Conjugate should negate bivectors
        conj_bivector = bivector.conjugate()
        self.assertEqual(conj_bivector.get(3), -1)
        self.assertEqual(conj_bivector.get(5), -2)
        self.assertEqual(conj_bivector.get(6), -3)
    
    def test_inverse(self):
        """Test the inverse operation."""
        # Create a scalar
        scalar = self.mk(e=2)
        # Inverse of scalar 2 should be 1/2
        inv_scalar = scalar.inv()
        self.assertEqual(inv_scalar.e, 0.5)
        
        # Create a vector
        vector = self.mk(**{1: 1})
        # Test that a vector can be inverted
        inv_vector = vector.inv()
        self.assertIsNotNone(inv_vector)
        
        # Test that vector * inv(vector) = 1
        prod = vector * inv_vector
        self.assertEqual(prod.grade(0).e, 1)
    
    def test_numpy_interop(self):
        """Test NumPy interoperability."""
        # Create an array of multivectors
        mvs = np.array([
            MultiVector.fromkeysvalues(self.vga3d, (1, 2), [i, i*2]) 
            for i in range(1, 4)
        ], dtype=object)
        
        # Test that the array has the right shape
        self.assertEqual(mvs.shape, (3,))
        
        # Test accessing elements
        self.assertEqual(mvs[0].get(1), 1)
        self.assertEqual(mvs[1].get(1), 2)
        self.assertEqual(mvs[2].get(1), 3)

if __name__ == "__main__":
    unittest.main()