from pathlib import Path
from setuptools import setup, find_packages

setup(
    name="kingdon",
    version="1.3.1",
    author="Your Name",
    author_email="your.email@example.com",
    description="Pythonic Geometric Algebra Package",
    long_description=(Path(__file__).parent / 'README.md').read_text(encoding='utf-8') if (Path(__file__).parent / 'README.md').exists() else '',
    long_description_content_type="text/markdown",
    url="https://github.com/tracyphasespace/kingdon",

    packages=['kingdon'],
    package_dir={'kingdon': 'src'},

    # === Add/modify this section for package data ===
    include_package_data=True, # Tells setuptools to look at MANIFEST.in and SCM data
    package_data={
        # If 'kingdon' is your package name and its modules are in 'SRC',
        # and 'graph.js' is also in 'SRC' alongside your .py files for that package:
        'kingdon': ['*.js', '*.css'], # Include all .js and .css files from the 'kingdon' package's root (which is SRC)
    },
    # ===============================================

    install_requires=[
        "numpy",
        "sympy",
        "anywidget" # <<< Make sure anywidget is listed as a dependency!
                    # Based on your pip list -v from turn 29, it was installed (anywidget 0.9.18)
                    # but it's good practice to have it in install_requires.
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Mathematics",
    ],
    python_requires='>=3.8',
    # Consider adding zip_safe=False if you have issues with package_data in editable mode,
    # though with modern setuptools and anywidget, it might not be necessary.
    # zip_safe=False,
)