#!/usr/bin/python3
# -*- coding: utf-8 -*-

""" A collection of utility functions related to the processing of SBS instrument data.

Functions:

    close_enough (np.ndarray, np.ndarray, int, float) -> bool
    plot (np.ndarray)
    percent_match (np.ndarray, np.ndarray) -> str

"""

# Native imports

# Third-party imports
import matplotlib.pyplot as plt
import numpy as np

# Sea-Bird imports

# Internal imports


def close_enough(
    test_values: np.ndarray,
    expected_values: np.ndarray,
    rounding_order: int,
    absolute_tolerance: float,
) -> bool:
    """Compares ndarrays, accounting for rounding errors introduced by SeaSoft.

    Frequently a float will be accurate to the 14th decimal, for example, but the
    15th decimal may toggle it over or under 5, causing it to round differently.
    This function adds or subtracts half of the least significant digit, then
    compares values to a given tolerance.

    Args:
        test_values (np.ndarray): values derived by this library
        expected_values (np.ndarray): values derived by SeaSoft
        rounding_order (int): The order that SeaSoft values were rounded to
        absolute_tolerance (float): values must be at least this close after
        adding or subtracting the rounding error

    Returns:
        bool: pass or fail aggregate
    """

    results = np.full(len(test_values), False)
    for n in range(len(test_values)):
        results[n] = (
            np.round(test_values[n], rounding_order) == expected_values[n]
            or np.isclose(
                test_values[n],
                expected_values[n] - 0.5 * 10**-rounding_order,
                rtol=0,
                atol=absolute_tolerance,
            )
            or np.isclose(
                test_values[n],
                expected_values[n] + 0.5 * 10**-rounding_order,
                rtol=0,
                atol=absolute_tolerance,
            )
        )
    return np.all(results)


def plot(**kwargs: np.ndarray):
    """Plots a dictionary of ndarrays

    Args:
        **kwargs (np.ndarray): the dictionary to plot
    """

    fig, ax = plt.subplots(figsize=(20, 10))
    for key in kwargs.keys():
        x = range(len(kwargs[key]))
        ax.plot(x, kwargs[key], label=key)
        ax.legend()
    plt.show()


def percent_match(x1: np.ndarray, x2: np.ndarray) -> str:
    """Calculates the extent that two arrays match.

    Args:
        x1 (np.ndarray): first array to be compared
        x2 (np.ndarray): second array to be compared

    Returns:
        str: a message with the percentage of matching elements
    """

    return f"{100 - (x1 != x2).sum()*100 / len(x1):0.2f}% match"
