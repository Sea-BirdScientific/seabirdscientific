"""EOS80 functions to support legacy conversions"""

import math
import warnings

import numpy as np
import seawater as sw
from scipy import stats

import seabirdscientific.constants as const


def buoyancy_eos80(
    temperature: np.ndarray,
    salinity: np.ndarray,
    pressure: np.ndarray,
    latitude: np.ndarray,
    longitude: np.ndarray,
    window_size: float,
    flag_value=const.FLAG_VALUE,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Calculates the 4 buoyancy values using the EOS-80 formula.

    Same as buoyancy, but uses eos80_conversion.bouyancy_frequency (the SBE
    Data Processing EOS-80 calculation) instead of the TEOS-10 formula.

    Data is expected to have already been binned via Bin_Average using
    decibar pressure bins. All arrays are expected to be the same length,
    except for latitude and longitude, which can be length 1.

    :param temperature: Temperature in ITS-90 degrees C
    :param salinity: Practical salinity in PSS-78 PSU
    :param pressure: Pressure in dbar
    :param latitude: latitude values. If length 1, gets applied to all values.
    :param longitude: longitude values. If length 1, gets applied to all values.
    :param window_size: window size to use. If this number is smaller than the
        binned window size, round up to a minimum of 3 scans.
    :param flag_value: Bad Flag value to use for marking bad scans.
        Defaults to -9.99e-29

    :return: a tuple of ndarrays including: buoyancy frequency squared,
        buoyancy frequency, stability, and scaled stability
    """

    _salinity, _temperature, _pressure, _latitude, _longitude = np.broadcast_arrays(
        salinity, temperature, pressure, latitude, longitude
    )

    # Get the original bin size using the second and third bin so we don't
    # have to worry about the surface bin
    original_bin_size = abs(_pressure[2] - _pressure[1])

    # Number of scans on either side of the median point, minimum 1
    scans_per_side = max(math.floor(window_size / original_bin_size / 2), 1)

    # create our result np.ndarrays with the flag value as default
    buoyancy_freq_squared = np.full(len(_temperature), flag_value)
    buoyancy_freq = np.full(len(_temperature), flag_value)
    stability = np.full(len(_temperature), flag_value)
    scaled_stability = np.full(len(_temperature), flag_value)

    for i in range(scans_per_side, len(_temperature) - scans_per_side):
        min_index = i - scans_per_side
        max_index = i + scans_per_side + 1  # + 1 because slicing excludes the max

        pressure_subset = _pressure[min_index:max_index]
        temperature_its_subset = _temperature[min_index:max_index]
        salinity_subset = _salinity[min_index:max_index]

        mean_pressure = np.mean(pressure_subset)
        # depth is negative below the surface (0 at the surface)
        depth = -sw.eos80.dpth(mean_pressure, _latitude[i])
        gravity = sw.eos80.g(_latitude[i], depth)

        n2 = bouyancy_frequency(temperature_its_subset, salinity_subset, pressure_subset, gravity)

        buoyancy_freq_squared[i] = n2
        if n2 >= 0:
            buoyancy_freq[i] = math.sqrt(n2) * 3600 / (2 * np.pi)
        else:
            # negative root of the absolute buoyancy squared to match seasoft
            buoyancy_freq[i] = -math.sqrt(abs(n2)) * 3600 / (2 * np.pi)
        stability[i] = n2 / gravity
        scaled_stability[i] = stability[i] * 1e8

    return (buoyancy_freq_squared, buoyancy_freq, stability, scaled_stability)


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
