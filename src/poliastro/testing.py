"""Testing utilities.

"""
import os.path

import pytest


def test(args=[]):
    """Initiate poliastro testing

    """
    pytest.main([os.path.dirname(os.path.abspath(__file__))] + args)
