from contextlib import contextmanager

from poliastro import jit


@contextmanager
def _fake_numba_import():
    # Black magic, beware
    # https://stackoverflow.com/a/2484402/554319
    import sys

    class FakeImportFailure:
        def __init__(self, modules):
            self.modules = modules

        def find_module(self, fullname, *args, **kwargs):
            if fullname in self.modules:
                raise ImportError('Debug import failure for %s' % fullname)

    fail_loader = FakeImportFailure(['numba'])

    import poliastro.jit
    from poliastro import jit
    del poliastro.jit
    del jit
    del sys.modules['poliastro.jit']

    del sys.modules['numba']

    sys.meta_path.insert(0, fail_loader)

    yield

    sys.meta_path.remove(fail_loader)
    import numba
    from poliastro import jit


def test_ijit_returns_same_function_without_args():
    def expected_foo():
        return True
    foo = jit.ijit(expected_foo)
    assert foo is expected_foo


def test_ijit_returns_same_function_with_args():
    def expected_foo():
        return True
    foo = jit.ijit(1)(expected_foo)
    assert foo is expected_foo


def test_no_numba_emits_warning(recwarn):
    with _fake_numba_import():
        from poliastro import jit

        assert len(recwarn) == 1
        w = recwarn.pop(UserWarning)
        assert issubclass(w.category, UserWarning)
        assert "Could not import numba package" in str(w.message)
