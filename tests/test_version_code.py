"""Validate the version code of the package."""

from validate_version_code import validate_version_code
from map4.__version__ import __version__


def test_version():
    """Test the package version."""
    assert validate_version_code(__version__)
