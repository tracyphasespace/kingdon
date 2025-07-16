#!/usr/bin/env python3
"""
Performance optimizations for geometric algebra operations.
Implements vectorization, multiprocessing, and GPU acceleration.
"""

import numpy as np
import multiprocessing as mp
from typing import List, Union, Optional
from functools import partial
import time

# Optional GPU support
try:
    import cupy as cp
    GPU_AVAILABLE = True
except ImportError:
    GPU_AVAILABLE = False
    cp = None

# Optional JIT compilation
try:
    import numba
    from numba import jit, cuda
    JIT_AVAILABLE = True
except ImportError:
    JIT_AVAILABLE = False
    numba = None

class PerformanceOptimizer:
    """
    Performance optimization utilities for geometric algebra operations.
    """
    
    def __init__(self, use_gpu=None, use_jit=None, num_processes=None):
        """
        Initialize performance optimizer.
        
        Args:
            use_gpu: Enable GPU acceleration (auto-detect if None)
            use_jit: Enable JIT compilation (auto-detect if None)
            num_processes: Number of processes for multiprocessing (auto-detect if None)
        """
        self.use_gpu = use_gpu if use_gpu is not None else GPU_AVAILABLE
        self.use_jit = use_jit if use_jit is not None else JIT_AVAILABLE
        self.num_processes = num_processes if num_processes is not None else mp.cpu_count()
        
        print(f"🚀 Performance Optimizer initialized:")
        print(f"   GPU acceleration: {'✅' if self.use_gpu else '❌'}")
        print(f"   JIT compilation: {'✅' if self.use_jit else '❌'}")
        print(f"   Multiprocessing: {self.num_processes} cores")
    
    def vectorize_coefficients(self, multivectors: List['MultiVector']) -> np.ndarray:
        """
        Convert list of multivectors to vectorized coefficient array.
        
        Args:
            multivectors: List of MultiVector objects
            
        Returns:
            NumPy array of shape (n_multivectors, max_coefficients)
        """
        if not multivectors:
            return np.array([])
        
        # Find maximum number of coefficients
        max_coeffs = max(len(mv._values) for mv in multivectors)
        
        # Create coefficient matrix
        coeffs = np.zeros((len(multivectors), max_coeffs))
        
        for i, mv in enumerate(multivectors):
            coeffs[i, :len(mv._values)] = mv._values
        
        return coeffs
    
    def batch_add(self, multivectors1: List['MultiVector'], 
                  multivectors2: List['MultiVector']) -> List['MultiVector']:
        """
        Vectorized addition of multivector arrays.
        
        Args:
            multivectors1: First array of multivectors
            multivectors2: Second array of multivectors
            
        Returns:
            List of sum multivectors
        """
        if len(multivectors1) != len(multivectors2):
            raise ValueError("Arrays must have same length")
        
        if not multivectors1:
            return []
        
        # Get algebra from first multivector
        algebra = multivectors1[0].algebra
        
        # Vectorize coefficients
        coeffs1 = self.vectorize_coefficients(multivectors1)
        coeffs2 = self.vectorize_coefficients(multivectors2)
        
        # Use GPU if available
        if self.use_gpu and coeffs1.size > 1000:  # Only use GPU for larger arrays
            coeffs1_gpu = cp.asarray(coeffs1)
            coeffs2_gpu = cp.asarray(coeffs2)
            result_coeffs_gpu = coeffs1_gpu + coeffs2_gpu
            result_coeffs = cp.asnumpy(result_coeffs_gpu)
        else:
            result_coeffs = coeffs1 + coeffs2
        
        # Convert back to multivectors
        results = []
        for i in range(len(multivectors1)):
            # Use the key structure from the first multivector (assuming similar structure)
            keys = multivectors1[i]._keys
            values = result_coeffs[i, :len(keys)]
            mv = algebra.multivector.__class__.fromkeysvalues(algebra, keys, values)
            results.append(mv)
        
        return results
    
    def batch_scalar_multiply(self, multivectors: List['MultiVector'], 
                             scalars: Union[float, List[float]]) -> List['MultiVector']:
        """
        Vectorized scalar multiplication of multivector arrays.
        
        Args:
            multivectors: Array of multivectors
            scalars: Scalar value or array of scalars
            
        Returns:
            List of scaled multivectors
        """
        if not multivectors:
            return []
        
        algebra = multivectors[0].algebra
        coeffs = self.vectorize_coefficients(multivectors)
        
        # Handle scalar broadcasting
        if isinstance(scalars, (int, float)):
            scalar_array = np.full(len(multivectors), scalars)
        else:
            scalar_array = np.array(scalars)
        
        # Use GPU if available
        if self.use_gpu and coeffs.size > 1000:
            coeffs_gpu = cp.asarray(coeffs)
            scalar_gpu = cp.asarray(scalar_array)
            result_coeffs_gpu = coeffs_gpu * scalar_gpu[:, np.newaxis]
            result_coeffs = cp.asnumpy(result_coeffs_gpu)
        else:
            result_coeffs = coeffs * scalar_array[:, np.newaxis]
        
        # Convert back to multivectors
        results = []
        for i in range(len(multivectors)):
            keys = multivectors[i]._keys
            values = result_coeffs[i, :len(keys)]
            mv = algebra.multivector.__class__.fromkeysvalues(algebra, keys, values)
            results.append(mv)
        
        return results
    
    def parallel_geometric_product(self, multivectors1: List['MultiVector'], 
                                  multivectors2: List['MultiVector']) -> List['MultiVector']:
        """
        Parallel computation of geometric products.
        
        Args:
            multivectors1: First array of multivectors
            multivectors2: Second array of multivectors
            
        Returns:
            List of geometric product results
        """
        if len(multivectors1) != len(multivectors2):
            raise ValueError("Arrays must have same length")
        
        if not multivectors1:
            return []
        
        # For small arrays, use sequential processing
        if len(multivectors1) < 100:
            return [mv1.algebra.gp(mv1, mv2) for mv1, mv2 in zip(multivectors1, multivectors2)]
        
        # Use multiprocessing for larger arrays
        def compute_gp(args):
            mv1, mv2 = args
            return mv1.algebra.gp(mv1, mv2)
        
        with mp.Pool(processes=self.num_processes) as pool:
            tasks = list(zip(multivectors1, multivectors2))
            results = pool.map(compute_gp, tasks)
        
        return results
    
    def benchmark_optimization(self, algebra, array_size=1000, num_trials=3):
        """
        Benchmark optimization techniques.
        
        Args:
            algebra: Algebra instance
            array_size: Size of test arrays
            num_trials: Number of benchmark trials
            
        Returns:
            Dictionary of benchmark results
        """
        print(f"\n🏁 BENCHMARKING OPTIMIZATIONS (array size: {array_size:,})")
        print("=" * 60)
        
        # Create test data
        np.random.seed(42)
        mvs1 = [algebra.vector(np.random.random(3)) for _ in range(array_size)]
        mvs2 = [algebra.vector(np.random.random(3)) for _ in range(array_size)]
        scalars = np.random.random(array_size)
        
        results = {}
        
        # Benchmark addition
        print("Testing Addition:")
        
        # Sequential addition
        times = []
        for _ in range(num_trials):
            start = time.perf_counter()
            seq_results = [mv1 + mv2 for mv1, mv2 in zip(mvs1, mvs2)]
            times.append(time.perf_counter() - start)
        seq_time = np.mean(times)
        
        # Vectorized addition
        times = []
        for _ in range(num_trials):
            start = time.perf_counter()
            vec_results = self.batch_add(mvs1, mvs2)
            times.append(time.perf_counter() - start)
        vec_time = np.mean(times)
        
        add_speedup = seq_time / vec_time if vec_time > 0 else float('inf')
        print(f"  Sequential: {seq_time:.4f}s ({array_size/seq_time:.0f} ops/sec)")
        print(f"  Vectorized: {vec_time:.4f}s ({array_size/vec_time:.0f} ops/sec)")
        print(f"  Speedup: {add_speedup:.1f}x")
        
        results['addition'] = {
            'sequential_time': seq_time,
            'vectorized_time': vec_time,
            'speedup': add_speedup
        }
        
        # Benchmark scalar multiplication
        print("\nTesting Scalar Multiplication:")
        
        # Sequential scalar multiplication
        times = []
        for _ in range(num_trials):
            start = time.perf_counter()
            seq_results = [mv * scalar for mv, scalar in zip(mvs1, scalars)]
            times.append(time.perf_counter() - start)
        seq_time = np.mean(times)
        
        # Vectorized scalar multiplication
        times = []
        for _ in range(num_trials):
            start = time.perf_counter()
            vec_results = self.batch_scalar_multiply(mvs1, scalars)
            times.append(time.perf_counter() - start)
        vec_time = np.mean(times)
        
        mul_speedup = seq_time / vec_time if vec_time > 0 else float('inf')
        print(f"  Sequential: {seq_time:.4f}s ({array_size/seq_time:.0f} ops/sec)")
        print(f"  Vectorized: {vec_time:.4f}s ({array_size/vec_time:.0f} ops/sec)")
        print(f"  Speedup: {mul_speedup:.1f}x")
        
        results['scalar_multiplication'] = {
            'sequential_time': seq_time,
            'vectorized_time': vec_time,
            'speedup': mul_speedup
        }
        
        # Benchmark geometric product (if array is not too large)
        if array_size <= 1000:
            print("\nTesting Geometric Product:")
            
            # Sequential geometric product
            times = []
            for _ in range(num_trials):
                start = time.perf_counter()
                seq_results = [algebra.gp(mv1, mv2) for mv1, mv2 in zip(mvs1[:100], mvs2[:100])]
                times.append(time.perf_counter() - start)
            seq_time = np.mean(times)
            
            # Parallel geometric product
            times = []
            for _ in range(num_trials):
                start = time.perf_counter()
                par_results = self.parallel_geometric_product(mvs1[:100], mvs2[:100])
                times.append(time.perf_counter() - start)
            par_time = np.mean(times)
            
            gp_speedup = seq_time / par_time if par_time > 0 else float('inf')
            print(f"  Sequential: {seq_time:.4f}s ({100/seq_time:.0f} ops/sec)")
            print(f"  Parallel: {par_time:.4f}s ({100/par_time:.0f} ops/sec)")
            print(f"  Speedup: {gp_speedup:.1f}x")
            
            results['geometric_product'] = {
                'sequential_time': seq_time,
                'parallel_time': par_time,
                'speedup': gp_speedup
            }
        
        return results

# JIT-compiled functions (if Numba is available)
if JIT_AVAILABLE:
    @jit(nopython=True)
    def jit_vector_add(coeffs1, coeffs2):
        """JIT-compiled vector addition."""
        return coeffs1 + coeffs2
    
    @jit(nopython=True)
    def jit_scalar_multiply(coeffs, scalars):
        """JIT-compiled scalar multiplication."""
        return coeffs * scalars[:, np.newaxis]

# GPU kernels (if CuPy is available)
if GPU_AVAILABLE:
    def gpu_batch_operations():
        """Example GPU operations using CuPy."""
        pass

def create_performance_enhanced_algebra(algebra_class):
    """
    Create a performance-enhanced version of an algebra class.
    
    Args:
        algebra_class: Original algebra class
        
    Returns:
        Enhanced algebra class with performance optimizations
    """
    class PerformanceEnhancedAlgebra(algebra_class):
        """Algebra class with performance optimizations."""
        
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.optimizer = PerformanceOptimizer()
        
        def batch_add(self, multivectors1, multivectors2):
            """Vectorized batch addition."""
            return self.optimizer.batch_add(multivectors1, multivectors2)
        
        def batch_scalar_multiply(self, multivectors, scalars):
            """Vectorized batch scalar multiplication."""
            return self.optimizer.batch_scalar_multiply(multivectors, scalars)
        
        def parallel_gp(self, multivectors1, multivectors2):
            """Parallel geometric product computation."""
            return self.optimizer.parallel_geometric_product(multivectors1, multivectors2)
        
        def benchmark_performance(self, array_size=1000):
            """Benchmark performance optimizations."""
            return self.optimizer.benchmark_optimization(self, array_size)
    
    return PerformanceEnhancedAlgebra

def main():
    """Test performance optimizations."""
    from src.algebra import Algebra
    
    print("🚀 TESTING PERFORMANCE OPTIMIZATIONS")
    print("=" * 60)
    
    # Create enhanced algebra
    EnhancedAlgebra = create_performance_enhanced_algebra(Algebra)
    ga = EnhancedAlgebra(p=3, q=0, r=0)
    
    # Run benchmarks
    results = ga.benchmark_performance(array_size=1000)
    
    print(f"\n🎯 OPTIMIZATION SUMMARY")
    print("=" * 60)
    for operation, data in results.items():
        print(f"{operation.title()}: {data['speedup']:.1f}x speedup")
    
    return results

if __name__ == "__main__":
    main()