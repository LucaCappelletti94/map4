#!/usr/bin/env python

import os
import re
from setuptools import find_packages, setup

here = os.path.abspath(os.path.dirname(__file__))

# Get the long description from the relevant file
with open(os.path.join(here, "README.md"), encoding="utf8") as f:
    long_description = f.read()


def find_version(*file_paths):
    """Read the version number from a source file."""
    with open(os.path.join(here, *file_paths), "r", encoding="utf8") as fp:
        version_file = fp.read()
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")


__version__ = find_version("map4", "__version__.py")

test_deps = [
    "pytest",
    "pytest-cov",
    "pytest-readme",
    "pandas",
    "validate_version_code",
]

extras = {
    "test": test_deps,
}

setup(
    name="MAP4",
    version=__version__,
    description="MinHashed AtomPair Fingerprint of Radius 2",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Alice Capecchi",
    author_email="alice.capecchi@outlook.it",
    license="MIT",
    url="https://github.com/reymond-group/map4",
    packages=find_packages(exclude=["contrib", "docs", "tests*"]),
    install_requires=[
        "mhfp",
        "numpy",
        "rdkit",
        "tqdm",
        "scikit-learn",
        "matplotlib",
    ],
    tests_require=test_deps,
    extras_require=extras,
    classifiers=[
        "Development Status :: 4 - Alpha",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
    ],
    entry_points={
        'console_scripts': [
            'map4 = map4.cli:main',  # CLI command and entry point
        ],
    },
)
