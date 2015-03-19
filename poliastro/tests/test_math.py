# coding: utf-8
import math as pymath

import pytest

from poliastro import math


def test_factorial():
    assert pymath.factorial(1) == math.factorial(1)
    assert pymath.factorial(2) == math.factorial(2)
    assert pymath.factorial(3) == math.factorial(3)
    assert pymath.factorial(5) == math.factorial(5)
    assert pymath.factorial(10) == math.factorial(10)
    assert pymath.factorial(20) == math.factorial(20)


def test_factorial_raises_ValueError_on_negative():
    with pytest.raises(ValueError):
        math.factorial(-1)


def test_factorial_raises_ValueError_on_non_integer():
    with pytest.raises(ValueError):
        math.factorial(1.5)


def test_factorial_raises_ValueError_on_overflow():
    with pytest.raises(ValueError):
        math.factorial(21)
