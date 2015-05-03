# coding: utf-8
"""Testing utilities.

"""
import os.path
import pytest


def test(coverage=False):
    args = [os.path.dirname(os.path.abspath(__file__))]
    if coverage:
        # Monkey patch jit to prevent numba.jit from shadowing coverage
        # figures (slow)
        import poliastro.jit
        poliastro.jit.jit = poliastro.jit.ijit

        args += ["--cov", "poliastro"]

    pytest.main(args)
