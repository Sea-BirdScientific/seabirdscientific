"""EOS80 functions to support legacy conversions"""

import warnings

import numpy as np
import seawater as sw
from scipy import stats


def bouyancy_frequency(
    temperature: np.ndarray,
    salinity: np.ndarray,
    pressure: np.ndarray,
    gravity: float,
):
    """Calculates an N^2 value (buoyancy frequency) for the given window
    of temperature, salinity, and pressure, at the given latitude.

    Expects temperature as ITS-90 temperature, salinity as practical
    salinity, and pressure as dbar, all of the same length. Performs the
    calculation following the SBE Data Processing formula using E0S-80
    calculations for potential temp and density

    :param temperature: ITS-90 temperature values for the given window
    :param salinity: PSS-78 practical salinity values for the given window
    :param pressure: pressure values for the given window
    :param gravity: gravity value

    :return: A single N^2 [Brunt-Väisälä (buoyancy) frequency]
    """

    db_to_pa = 1e4

    # Wrap these as a length-1 array
    pressure_bar = np.array([np.mean(pressure)])
    temperature_bar = np.array([np.mean(temperature)])
    salinity_bar = np.array([np.mean(salinity)])

    # EOS-80 density from the seawater library. seawater.dens returns full
    # density; subtract 1000 to keep the sigma convention the SBE Data
    # Processing formula (and the SeaSoft reference) are built around.
    rho_bar = np.atleast_1d(sw.dens(salinity_bar, temperature_bar, pressure_bar))[0] - 1000.0

    # EOS-80 potential temperature and density, referenced to the window mean.
    theta = sw.ptmp(salinity, temperature, pressure, pressure_bar[0])
    v_vals = 1.0 / (np.atleast_1d(sw.dens(salinity, theta, pressure_bar[0])) - 1000.0)

    # Estimate vertical gradient of specific volume
    dvdp_result = stats.linregress(pressure, v_vals)

    # Compute EOS-80 N2 combining computed average density and vertical gradient
    # we index into v_bar, alpha_bar, and beta_bar as they are all arrays of len 1
    n2 = 0 - (rho_bar**2 * gravity**2 * dvdp_result.slope / db_to_pa)
    return n2


def density(
    salinity: np.ndarray,
    temperature: np.ndarray,
    pressure: np.ndarray,
) -> np.ndarray:
    """EOS-80 density calculation.

    Delegates to the seawater library (seawater.dens), returning sigma
    (density - 1000) to match the SBE Data Processing convention.

    :param salinity: salinity data
    :param temperature: temperature data
    :param pressure: pressure data

    :return: resulting density data
    """

    warnings.warn(
        "eos80_conversion.density is deprecated; use the seawater library "
        "(seawater.dens) instead.",
        DeprecationWarning,
        stacklevel=2,
    )

    # seawater.dens returns full density; subtract 1000 to keep the sigma
    # convention used by SBE Data Processing (and the SeaSoft reference).
    return np.atleast_1d(sw.dens(salinity, temperature, pressure)) - 1000.0


def potential_temperature(
    salinity: np.ndarray,
    temperature: np.ndarray,
    pressure: np.ndarray,
    mean_pressure: np.ndarray,
) -> np.ndarray:
    """EOS-80 potential temperature calculation.

    Delegates to the seawater library (seawater.ptmp).

    :param salinity: sainity data
    :param temperature: temperature data
    :param pressure: subset pressure data
    :param mean_pressure: pressure data

    :return: calculated potential temperature data
    """

    warnings.warn(
        "eos80_conversion.potential_temperature is deprecated; use the seawater "
        "library (seawater.ptmp) instead.",
        DeprecationWarning,
        stacklevel=2,
    )

    return np.atleast_1d(sw.ptmp(salinity, temperature, pressure, mean_pressure))


def adiabatic_temperature_gradient(
    salinity: np.ndarray,
    temperature: np.ndarray,
    pressure: np.ndarray,
) -> np.ndarray:
    """EOS-80 adiabatic lapse rate calculation.

    This was ported from CSharedCalc::ATG()

    :param salinity: salinity data
    :param temperature: temperature data
    :param pressure: pressure data

    :return: the resulting adiabatic lapse rate
    """

    s, t, p = np.broadcast_arrays(salinity, temperature, pressure)

    ds = s - 35.0
    atg = (
        (
            ((-2.1687e-16 * t + 1.8676e-14) * t - 4.6206e-13) * p
            + (
                (2.7759e-12 * t - 1.1351e-10) * ds
                + ((-5.4481e-14 * t + 8.733e-12) * t - 6.7795e-10) * t
                + 1.8741e-8
            )
        )
        * p
        + (-4.2393e-8 * t + 1.8932e-6) * ds
        + ((6.6228e-10 * t - 6.836e-8) * t + 8.5258e-6) * t
        + 3.5803e-5
    )

    return atg
