"""Utility functions related to processing SBS instrument data."""

# Native imports

# Third-party imports
import matplotlib.pyplot as plt
import numpy as np
from line_profiler import LineProfiler

# Sea-Bird imports

# Internal imports
from seabirdscientific.processing import FLAG_VALUE


def close_enough(
    test_values: np.ndarray,
    expected_values: np.ndarray,
    rounding_order: int,
    absolute_tolerance: float,
) -> bool:
    """Compares ndarrays, ignoring differences due to rounding or
    truncating the least significant digit. This is only for comparing
    data to legacy software

    Occasionally a float will be accurate to some decimal, but the next
    lower decimal may toggle it over or under 5, causing it to round
    differently if it's rounded or truncated. This function adds and
    subtracts half of the least significant digit, then compares values
    of a given tolerance to either result

    :param test_values: values derived by this library
    :param expected_values: values derived by SeaSoft
    :param rounding_order: The order that CNV values were rounded to
    :param absolute_tolerance: values must be at least this close after
        adding or subtracting the rounding error

    :return: pass or fail aggregate
    """

    results = np.full(len(test_values), False)
    for n, v in enumerate(test_values):
        results[n] = (
            np.round(v, rounding_order) == expected_values[n]
            or np.isclose(
                v,
                expected_values[n] - 0.5 * 10**-rounding_order,
                rtol=0,
                atol=absolute_tolerance,
            )
            or np.isclose(
                v,
                expected_values[n] + 0.5 * 10**-rounding_order,
                rtol=0,
                atol=absolute_tolerance,
            )
        )
    return bool(np.all(results))


def plot(**kwargs: np.ndarray):
    """Plots a dictionary of ndarrays

    :param kwargs: the dictionary to plot
    """

    _, ax = plt.subplots(figsize=(20, 10))
    for key, value in kwargs.items():
        x = range(len(value))
        ax.plot(x, value, label=key)
        ax.legend()
    plt.show()


def percent_match(x1: np.ndarray, x2: np.ndarray) -> str:
    """Calculates the extent that two arrays match.

    :param x1: first array to be compared
    :param x2: second array to be compared

    :return: a message with the percentage of matching elements
    """

    return f"{100 - (x1 != x2).sum() * 100 / len(x1):0.2f}% match"


def get_tolerance(data: np.ndarray):
    """Checks the first 10 values of an array, gets the longest decimal
    length to the right of the decimal, and returns 1/10^length. Used
    for unit tests so results can be compared with a variable tolerance

    :param data: a numpy array of numbers
    :return: the number of significant digits to the right of the decimal
    """
    decimal_lengths = [0]
    for n in range(len(data)):
        if not np.isnan(data[n]) and data[n] != FLAG_VALUE:
            decimal_lengths.append(len(f"{data[n]}".split(".")[1]))
            if len(decimal_lengths) >= 10:
                break
    return 1 / 10 ** max(decimal_lengths)


def profile(fun):
    """Decorator for profiling long running functions during
    development. Add @profile above the function to be measured, then
    call the function in a script to get a line by line report printed
    to the console. Remove the decorator when done

    :param fun: This is implicitly the function below the decorator
    """

    def wrapper(*args, **kwargs):
        lp = LineProfiler()
        lp.add_function(fun)
        lp.runctx("result = func(*args, **kwargs)", globals(), locals())
        lp.print_stats(output_unit=1e-6)
        return locals()["result"]

    return wrapper
