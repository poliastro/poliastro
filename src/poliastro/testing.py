# coding: utf-8
"""Testing utilities.

"""
import os.path
import pytest


def test():
    pytest.main(os.path.dirname(os.path.abspath(__file__)))
