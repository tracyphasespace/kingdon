"""
Graph Widget Module for Geometric Algebra Visualization
======================================================

This module provides visualization capabilities for Geometric Algebra objects using
a JavaScript-based renderer (typically ganja.js). It's designed to work in Jupyter
notebooks and other compatible environments.

The module contains:
1. Helper functions for encoding MultiVector objects into JSON-serializable format
2. The GraphWidget class that serves as the bridge between Python and JavaScript

Key components:
- `walker`: Recursively flattens nested generators and iterables
- `encode`: Converts MultiVector objects and other data structures into JSON-serializable format
- `GraphWidget`: Main widget class that handles rendering and interaction

Example:
    >>> from kingdon.algebra import Algebra
    >>> from kingdon.graph import GraphWidget
    >>> alg = Algebra(p=3, q=0, r=1)  # 3D projective geometric algebra
    >>> v1 = alg.multivector(values=[1, 0, 0, 0], grades=(1,))
    >>> v2 = alg.multivector(values=[0, 1, 0, 0], grades=(1,))
    >>> # Render two vectors
    >>> GraphWidget(algebra=alg, raw_subjects=[v1, v2], options={"grid": True})
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping
from functools import cached_property
from types import GeneratorType
import pathlib
from typing import Any, Dict, Generator, List, Optional, Sequence, Set, Tuple, Type, Union, cast

import anywidget
import traitlets
import numpy as np

from kingdon.multivector import MultiVector


# Define the tree types for recursive encoding
TREE_TYPES: Tuple[Type, ...] = (list, tuple)


def walker(encoded_generator: Any, tree_types: Tuple[Type, ...] = TREE_TYPES) -> List[Any]:
    """
    Recursively traverse and flatten nested iterables into a single list.
    
    This function handles generators, lists, tuples, and other iterables,
    converting them all into a flat list structure suitable for serialization.
    
    Args:
        encoded_generator: A generator or iterable with potentially nested structure
        tree_types: Types to be treated as trees and recursively processed
        
    Returns:
        A flat list containing all elements from the nested structure
        
    Example:
        >>> walker([1, 2, [3, 4, (5, 6)]])
        [1, 2, [3, 4, [5, 6]]]
    """
    result = []
    for item in encoded_generator:
        if isinstance(item, GeneratorType):
            result.extend(walker(item, tree_types))
        elif isinstance(item, tree_types):
            result.append(walker(item, tree_types))
        else:
            result.append(item)
    return result


def encode(obj: Any, tree_types: Tuple[Type, ...] = TREE_TYPES, root: bool = False) -> Generator[Any, None, None]:
    """
    Encode objects into JSON-serializable structures for visualization.
    
    This function handles:
    - MultiVector objects (converting their internal data)
    - Nested iterables (recursively processing their elements)
    - Callable objects (executing them and encoding their results)
    
    Args:
        obj: The object to encode (MultiVector, iterable, callable, etc.)
        tree_types: Types to be treated as trees and recursively processed
        root: If True and obj is a tree type, yields encoded elements directly
        
    Yields:
        Encoded objects suitable for JSON serialization
        
    Example:
        >>> mv = algebra.multivector(values=[1, 0, 0], grades=(1,))
        >>> list(encode(mv))
        [{'mv': [1, 0, 0], 'keys': (1, 2, 4)}]
        
        >>> list(encode([mv, 5]))
        [[{'mv': [1, 0, 0], 'keys': (1, 2, 4)}, 5]]
    """
    if root and isinstance(obj, tree_types):
        # If root and obj is a tree, yield encoded elements directly
        yield from (encode(value, tree_types) for value in obj)
    elif isinstance(obj, tree_types):
        # For tree types, recursively encode all elements and maintain structure
        yield obj.__class__(encode(value, tree_types) for value in obj)
    elif isinstance(obj, MultiVector) and hasattr(obj, 'shape') and len(getattr(obj, 'shape', ())) > 1:
        # For multidimensional multivectors, iterate over each individual multivector
        if hasattr(obj, 'itermv'):
            yield from (encode(value) for value in obj.itermv())
        else:
            # Fallback for older versions without itermv method
            for idx in np.ndindex(obj.shape):
                yield from encode(obj[idx])
    elif isinstance(obj, MultiVector):
        # For regular multivectors, encode values and optionally keys
        if isinstance(obj._values, np.ndarray):
            # For NumPy array values, convert to bytes for efficient transfer
            values = obj._values.tobytes()
        else:
            # For other value types, make a copy
            values = list(obj._values) if hasattr(obj._values, '__iter__') else obj._values
            
        # Include keys only if the multivector isn't full
        if len(obj) != len(obj.algebra):
            yield {'mv': values, 'keys': obj._keys}
        else:
            yield {'mv': values}
    elif isinstance(obj, Callable):
        # For callable objects, execute them and encode the result
        yield from encode(obj(), tree_types)
    else:
        # For other objects, yield them directly
        yield obj


class GraphWidget(anywidget.AnyWidget):
    """
    Widget for visualizing Geometric Algebra objects in Jupyter notebooks.
    
    This widget provides an interactive visualization of MultiVector objects,
    with support for dragging, rotation, and other interactive features.
    
    Attributes:
        algebra: The Geometric Algebra instance for the visualization
        options: Dictionary of rendering options passed to the JavaScript renderer
        raw_subjects: Original objects to visualize
        signature: The signature of the algebra as a list of integers
        cayley: The Cayley table for multiplication in the algebra
        key2idx: Mapping from binary keys to indices
        subjects: Encoded subjects ready for rendering
        draggable_points: MultiVector objects that can be dragged
        draggable_points_idxs: Indices of draggable points in subjects list
        
    Example:
        >>> # Create a visualization of a vector and a plane
        >>> vector = create_vector(algebra, [1, 0, 0])
        >>> plane = create_bivector(algebra, [0, 0, 1])
        >>> GraphWidget(
        ...     algebra=algebra,
        ...     raw_subjects=[vector, plane],
        ...     options={"grid": True, "scale": 1.5}
        ... )
    """

    # Path to the JavaScript module for rendering
    _esm = pathlib.Path(__file__).parent / "graph.js"

    # Core properties
    algebra = traitlets.Instance(klass=object, allow_none=False, 
                                help="The Geometric Algebra instance")
    options = traitlets.Dict(default_value={}, help="Rendering options").tag(sync=True)

    # Properties synchronized with JavaScript
    signature = traitlets.List([], help="The algebra signature").tag(sync=True)
    cayley = traitlets.List([], help="Multiplication table").tag(sync=True)
    key2idx = traitlets.Dict({}, help="Key to index mapping").tag(sync=True)

    # Visualization properties
    raw_subjects = traitlets.List([], help="Original objects to visualize")
    pre_subjects = traitlets.List([], help="Processed subjects for visualization")
    draggable_points = traitlets.List([], help="Points that can be dragged").tag(sync=True)
    draggable_points_idxs = traitlets.List([], help="Indices of draggable points").tag(sync=True)
    subjects = traitlets.List([], help="Encoded subjects for visualization").tag(sync=True)

    def __init__(self, algebra: Any, raw_subjects: List[Any], 
                 options: Optional[Dict[str, Any]] = None, **kwargs: Any):
        """
        Initialize the GraphWidget with algebra, subjects to visualize, and options.
        
        Args:
            algebra: The Geometric Algebra instance
            raw_subjects: List of objects to visualize (MultiVectors, callables, etc.)
            options: Dictionary of visualization options
            **kwargs: Additional arguments passed to the parent class
        """
        options = options or {}
        super().__init__(algebra=algebra, raw_subjects=raw_subjects, options=options, **kwargs)
        # Register message handler for JavaScript interactions
        self.on_msg(self._handle_custom_msg)

    def _handle_custom_msg(self, data: Dict[str, Any], buffers: Any) -> None:
        """
        Handle messages from the JavaScript frontend.
        
        This allows for bidirectional communication between Python and JavaScript.
        
        Args:
            data: Message data from the frontend
            buffers: Binary buffers attached to the message
        """
        if data.get("type") == "update_mvs":
            # Update the subjects when dragging occurs in the frontend
            self.subjects = self.get_subjects()

    def _get_pre_subjects(self) -> Any:
        """
        Process raw_subjects to prepare them for visualization.
        
        If there's a single callable subject (not a MultiVector), call it and use its result.
        
        Returns:
            Processed subjects ready for encoding
        """
        if len(self.raw_subjects) == 1:
            s = self.raw_subjects[0]
            if not isinstance(s, MultiVector) and isinstance(s, Callable):
                pre_subjects = s()
                if not isinstance(pre_subjects, TREE_TYPES):
                    pre_subjects = [pre_subjects]
                return pre_subjects
        return self.raw_subjects

    @traitlets.default('key2idx')
    def get_key2idx(self) -> Dict[Any, int]:
        """
        Create a mapping from binary keys to indices.
        
        Returns:
            Dictionary mapping keys to their positions
        """
        return {k: i for i, k in enumerate(self.algebra.canon2bin.values())}

    @traitlets.default('signature')
    def get_signature(self) -> List[int]:
        """
        Get the algebra signature as a list of integers.
        
        Returns:
            List of integers representing the signature
        """
        return [int(s) for s in self.algebra.signature]

    @traitlets.default('cayley')
    def get_cayley(self) -> List[List[str]]:
        """
        Generate the Cayley table for multiplication in the algebra.
        
        Returns:
            2D list representing the multiplication table
        """
        # Check if algebra has a cayley attribute
        if hasattr(self.algebra, 'cayley'):
            return [
                [(s if (s := self.algebra.cayley[eJ, eI])[-1] != 'e' else f"{s[:-1]}1")
                 for eI in self.algebra.canon2bin]
                for eJ in self.algebra.canon2bin
            ]
        
        # Fallback implementation if cayley isn't directly available
        d = self.algebra.d
        n = 2**d
        cayley = []
        for i in range(n):
            row = []
            for j in range(n):
                # Determine the sign and result of multiplying basis blades i and j
                sign = self.algebra.signs.get((i, j), 0)
                if sign == 0:
                    row.append("0")
                else:
                    key = i ^ j  # XOR for geometric product
                    prefix = "-" if sign < 0 else ""
                    suffix = "" if key > 0 else "1"
                    row.append(f"{prefix}e{key}{suffix}")
            cayley.append(row)
        return cayley

    @traitlets.default('pre_subjects')
    def get_pre_subjects(self) -> Any:
        """
        Get the preprocessed subjects.
        
        Returns:
            Processed subjects ready for encoding
        """
        return self._get_pre_subjects()

    @traitlets.default('subjects')
    def get_subjects(self) -> List[Any]:
        """
        Get the final encoded subjects for visualization.
        
        Returns:
            List of encoded subjects ready for JavaScript rendering
        """
        return walker(encode(self._get_pre_subjects(), root=True))

    @traitlets.default('draggable_points')
    def get_draggable_points(self) -> Any:
        """
        Identify draggable points from the subjects.
        
        In PGA (r=1, d=3 or d=4), only points of specific grades are draggable.
        
        Returns:
            List of encoded draggable points
        """
        d = self.algebra.d
        points = [s for s in self.pre_subjects if isinstance(s, MultiVector)]
        
        # In PGA, only points of specific grades are draggable
        if self.algebra.r == 1 and (d == 3 or d == 4):
            points = [p for p in points if p.grades == (d - 1,)]
            
        return walker(encode(points))

    @traitlets.default('draggable_points_idxs')
    def get_draggable_points_idxs(self) -> List[int]:
        """
        Get the indices of draggable points in the subjects list.
        
        Returns:
            List of indices for draggable points
        """
        d = self.algebra.d
        if self.algebra.r == 1 and (d == 3 or d == 4):
            return [j for j, s in enumerate(self.pre_subjects)
                    if isinstance(s, MultiVector) and s.grades == (d - 1,)]
        return [j for j, s in enumerate(self.pre_subjects) if isinstance(s, MultiVector)]

    @traitlets.observe('draggable_points')
    def _observe_draggable_points(self, change: Dict[str, Any]) -> None:
        """
        Update subjects when draggable points change.
        
        This is triggered when points are dragged in the frontend.
        
        Args:
            change: Dictionary containing the new values
        """
        self.inplacereplace(self.pre_subjects, list(zip(self.draggable_points_idxs, change['new'])))
        self.subjects = self.get_subjects()

    @traitlets.validate("options")
    def _valid_options(self, proposal: Dict[str, Any]) -> Dict[str, Any]:
        """
        Validate and process the options dictionary.
        
        Special handling for options like 'camera' that need encoding.
        
        Args:
            proposal: The proposed new value
            
        Returns:
            Processed options dictionary
        """
        options = proposal['value']
        if 'camera' in options:
            options['camera'] = list(encode(options['camera']))[0]
        return options

    def inplacereplace(self, old_subjects: List[Any], new_subjects: List[Tuple[int, Dict[str, Any]]]) -> None:
        """
        Update values in old_subjects based on changes in new_subjects.
        
        This performs an in-place update to avoid creating new objects unnecessarily.
        
        Args:
            old_subjects: Original list of subjects
            new_subjects: List of (index, new_value) pairs
        """
        for j, new_subject in new_subjects:
            old_subject = old_subjects[j]
            old_vals = old_subject._values
            new_vals = new_subject['mv']
            
            if len(old_vals) == len(self.algebra):
                # Full multivector: update all values
                for i, val in enumerate(new_vals):
                    if old_vals[i] != val:
                        old_vals[i] = val
            else:
                # Partial multivector: update values by key
                for i, k in enumerate(old_subject._keys):
                    val = new_vals[self.key2idx[k]]
                    if old_vals[i] != val:
                        old_vals[i] = val