# coding: utf-8
from poliastro import jit


def test_ijit_returns_same_function_without_args():
    expected_foo = lambda: True
    foo = jit.ijit(expected_foo)
    assert foo is expected_foo


def test_ijit_returns_same_function_with_args():
    expected_foo = lambda: True
    foo = jit.ijit(1)(expected_foo)
    assert foo is expected_foo
