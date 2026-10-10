"""Import dependency unit tests."""

import os
import subprocess
import sys
import textwrap

import pytest

# Plotting and profiling packages, which the data modules shouldn't need
PLOTTING_AND_PROFILING_MODULES = ["matplotlib", "plotly", "line_profiler"]

CORE_MODULES = [
    "seabirdscientific.cal_coefficients",
    "seabirdscientific.constants",
    "seabirdscientific.contour",
    "seabirdscientific.conversion",
    "seabirdscientific.eos80_conversion",
    "seabirdscientific.eos80_processing",
    "seabirdscientific.instrument_data",
    "seabirdscientific.interpret_sbs_variable",
    "seabirdscientific.processing",
    "seabirdscientific.utils",
]


def run_python(
    code: str, block_plotting_and_profiling: bool = False
) -> subprocess.CompletedProcess:
    """Runs code in a fresh interpreter

    :param code: python source to run
    :param block_plotting_and_profiling: if True, importing matplotlib,
        plotly or line_profiler raises ImportError, as if they weren't
        installed
    :return: the completed process
    """
    blocker = ""
    if block_plotting_and_profiling:
        blocker = "import sys\n" + "".join(
            f"sys.modules[{name!r}] = None\n" for name in PLOTTING_AND_PROFILING_MODULES
        )
    return subprocess.run(
        [sys.executable, "-c", blocker + textwrap.dedent(code)],
        capture_output=True,
        text=True,
        check=False,
        env={**os.environ, "MPLBACKEND": "Agg"},
    )


class TestImportDependencies:
    @pytest.mark.parametrize("module", CORE_MODULES)
    def test_core_module_imports_without_plotting_or_profiling(self, module):
        result = run_python(f"import {module}", block_plotting_and_profiling=True)
        assert result.returncode == 0, result.stderr

    def test_utils_plot_still_works(self):
        result = run_python(
            """
            import numpy as np
            from seabirdscientific.utils import plot
            plot(x=np.zeros(3))
            """
        )
        assert result.returncode == 0, result.stderr

    def test_utils_profile_still_decorates(self):
        # Only decorates: calling a profiled function hits a separate issue on Python 3.13+,
        # where the wrapper's locals()["result"] lookup fails
        result = run_python(
            """
            from seabirdscientific.utils import profile

            @profile
            def add(a, b):
                return a + b
            """
        )
        assert result.returncode == 0, result.stderr
