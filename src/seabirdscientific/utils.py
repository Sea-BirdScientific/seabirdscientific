"""Utility functions related to processing SBS instrument data."""

# Native imports
import warnings
from enum import EnumMeta

# Third-party imports
import matplotlib.pyplot as plt
import numpy as np
from line_profiler import LineProfiler

# Sea-Bird imports
import seabirdscientific.constants as const


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


def get_tolerance(data: np.ndarray, flag_value=-9.99e-29):
    """Checks the first 10 values of an array, gets the longest decimal
    length to the right of the decimal, and returns 1/10^length. Used
    for unit tests so results can be compared with a variable tolerance

    :param data: a numpy array of numbers
    :return: the number of significant digits to the right of the decimal
    """
    decimal_lengths = [0]
    for n in range(len(data)):
        if not np.isnan(data[n]) and data[n] != flag_value:
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
        lp.runctx("result = fun(*args, **kwargs)", globals(), locals())
        lp.print_stats(output_unit=1e-6)
        return locals()["result"]

    return wrapper


class WarnAllMembersMeta(EnumMeta):
    def __getattribute__(cls, name):
        obj = super().__getattribute__(name)
        if isinstance(obj, cls):
            warnings.warn(f"{cls.__name__}.{name} is deprecated", DeprecationWarning, stacklevel=2)
        return obj


def compute_scaled_temperature(temperature: np.ndarray) -> np.ndarray:
    return np.log((const.KELVIN_OFFSET_25C - temperature) / (const.KELVIN_OFFSET_0C + temperature))


def compute_ln_salinity_correction(temperature: np.ndarray, salinity: np.ndarray) -> np.ndarray:
    """Compute natural logarithm of the salinity correction for Garcia and Gordon
    Oxygen Solubility. Also applicable to SBE 63 Oxygen

    :param temperature: Temperature in degrees Celsius
    :param salinity: Salinity in PSU

    :return: Natural logarithm of the salinity correction"""
    sol_b0 = -6.24523e-3
    sol_b1 = -7.37614e-3
    sol_b2 = -1.0341e-2
    sol_b3 = -8.17083e-3
    sol_c0 = -4.88682e-7

    ts = compute_scaled_temperature(temperature)
    s_corr = (
        salinity * (sol_b0 + sol_b1 * ts + sol_b2 * ts**2 + sol_b3 * ts**3) + sol_c0 * salinity**2
    )
    return s_corr


def compute_rolling_average(
    compute_var: np.ndarray,
    window_size: float,
    sample_interval: float,
    modification_fn: callable = None,
) -> np.ndarray:
    """Computes a rolling average of the given variable over the specified window size.
    Averages with equal number values on either side of the center of each window.

    :param compute_var: The variable to compute the rolling average for.
    :param window_size: The size of the rolling window (in seconds).
    :param sample_interval: The time interval between samples (in seconds).
    :param modification_fn: Optional function to modify the computed rolling average.

    :return: An array containing the rolling average values.
    """
    if window_size <= 0:
        raise ValueError("Window size must be a positive integer.")

    # Calculate the number of samples in the rolling window
    num_samples = int(window_size / sample_interval)

    # Determine padding needed for both ends
    pad_before = num_samples // 2
    pad_after = num_samples - 1 - pad_before

    # Pad the array using the edge values
    # This prevents the ends from dropping off or pulling toward zero
    padded_data = np.pad(compute_var, (pad_before, pad_after), mode="edge")

    # Compute the rolling average using numpy's convolve function
    weights = np.ones(num_samples) / num_samples
    rolling_avg = np.convolve(padded_data, weights, mode="valid")

    # Apply the modification function if provided
    if modification_fn is not None:
        rolling_avg = modification_fn(rolling_avg)

    return rolling_avg
