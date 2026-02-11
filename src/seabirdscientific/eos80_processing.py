"""EOS80 functions to support legacy processing"""

import warnings
import numpy as np

import seabirdscientific.eos80_conversion as ec


def bouyancy_frequency(
    temp_ITS_subset: np.ndarray,
    salinity_prac_subset: np.ndarray,
    pressure_dbar_subset: np.ndarray,
    gravity: float,
):
    """Calculates an N^2 value (buoyancy frequency) for the given window
    of temperature, salinity, and pressure, at the given latitude.

    Expects temperature as ITS-90 temperature, salinity as practical
    salinity, and pressure as dbar, all of the same length. Performs the
    calculation following the SBE Data Processing formula using E0S-80
    calculations for potential temp and density

    :param temperature: ITS-90 temperature values for the given window
    :param salinity: practical salinity values for the given window
    :param pressure: pressure values for the given window
    :param gravity: gravity value

    :return: A single N^2 [Brunt-Väisälä (buoyancy) frequency]
    """

    warnings.warn("Deprecated. Use eos80_conversion.buoyancy_frequency", DeprecationWarning)

    return ec.bouyancy_frequency(
        temp_ITS_subset,
        salinity_prac_subset,
        pressure_dbar_subset,
        gravity,
    )


def density(
    s0: np.ndarray,
    t: np.ndarray,
    p0: np.ndarray,
) -> np.ndarray:
    """Deprecated. Use eos80_conversion.density"""

    warnings.warn("Deprecated. Use eos80_conversion.density", DeprecationWarning)

    return ec.density(s0, t, p0)


def potential_temperature(
    s: np.ndarray, t0: np.ndarray, p0: np.ndarray, pr: np.ndarray
) -> np.ndarray:
    """Deprecated. Use eos80_conversion.potential_temperature"""

    warnings.warn("Deprecated. Use eos80_conversion.potential_temperature", DeprecationWarning)

    return ec.potential_temperature(s, t0, p0, pr)


def adiabatic_temperature_gradient(
    s: np.ndarray,
    t: np.ndarray,
    p: np.ndarray,
) -> np.ndarray:
    """Deprecated. Use eos80_conversion.adiabatic_temperature_gradient"""

    warnings.warn("Deprecated. Use eos80_conversion.adiabatic_temperature_gradient", DeprecationWarning)

    return ec.adiabatic_temperature_gradient(s, t, p)
