import inspect

from astropy import units as u

from poliastro.twobody.decorators import state_from_vector


def fun(t, ss):
    return 1


def decorated_fun(t, u_, k):
    pass


def test_decorator_has_correct_signature():
    expected_signature = inspect.getfullargspec(decorated_fun)

    new_fun = state_from_vector(fun)

    assert inspect.getfullargspec(new_fun).args == expected_signature.args


def test_decorated_function_calls_rvstate():
    _t = 0
    _u = [1, 0, 0, 0, 1, 0]
    _k = 1

    new_fun = state_from_vector(fun)

    assert new_fun(_t, _u, _k) == 1
