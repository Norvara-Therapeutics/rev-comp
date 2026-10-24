"""The examples in the package documentation are executed here.

Running them through :mod:`doctest` against the imported package keeps them
correct wherever the package is installed from, rather than collecting the
source files by path, which cannot work when the import resolves elsewhere.
"""

import doctest
import importlib
import pkgutil

import revcomp


def submodules():
    """Every module in the package, so a new one is covered automatically."""
    yield revcomp
    for info in pkgutil.iter_modules(revcomp.__path__, revcomp.__name__ + "."):
        yield importlib.import_module(info.name)


def test_docstring_examples_run_and_pass():
    attempted = 0
    for module in submodules():
        failed, ran = doctest.testmod(module, verbose=False)
        assert failed == 0, (
            f"{module.__name__}: {failed} of {ran} docstring example(s) failed"
        )
        attempted += ran
    assert attempted, "no docstring examples were found"
