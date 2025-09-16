#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""A collection of utility functions related to the processing of SBS
instrument data.
"""
# Functions:
#   close_enough (np.ndarray, np.ndarray, int, float) -> bool
#   plot (np.ndarray)
#   percent_match (np.ndarray, np.ndarray) -> str

# Native imports

# Third-party imports
import matplotlib.pyplot as plt
import numpy as np

# Sea-Bird imports

# Internal imports
from .processing import FLAG_VALUE


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

    return f"{100 - (x1 != x2).sum()*100 / len(x1):0.2f}% match"


def get_decimal_length(data: np.ndarray):
    """Checks the first 10 values of an array and returns the longest
    decimal length to the right of the decimal. Used for unit tests so
    a variable tolerance can be used to compare results

    :param data: a numpy array of numbers
    :return: the number of significant digits to the right of the decimal
    """
    decimal_lengths = [0]
    for n in range(len(data)):
        if not np.isnan(data[n]) and data[n] != FLAG_VALUE:
            decimal_lengths.append(len(f"{data[n]}".split(".")[1]))
            if len(decimal_lengths) >= 10:
                break
    return max(decimal_lengths)
