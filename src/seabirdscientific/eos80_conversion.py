"""EOS80 functions to support legacy conversions"""

import math
import warnings

import numpy as np
import seawater as sw
from scipy import stats
from typing import Literal

import seawater as sw

from . import constants as const

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


def derive_potential_temperature_anomaly(
    salinity: np.ndarray,
    temperature: np.ndarray,
    pressure: np.ndarray,
    a0: float = 0.0,
    a1: float = 0.0,
    a1multiplier: Literal["salinity", "sigma-theta"] = "salinity",
) -> np.ndarray:
    """Derive potential temperature anomaly from practical salinity,
    temperature, and pressure using EOS-80 calculations.

    Calculates potential temperature at reference pressure (0 dbar) and
    applies anomaly corrections based on the provided coefficients a0 and a1.
    Uses EOS-80 formulations for potential temperature and density calculations.

    :param salinity: Practical salinity in PSU (ndarray)
    :param temperature: Temperature in ITS-90 degrees C (ndarray)
    :param pressure: Pressure in decibars (ndarray)
    :param a0: Anomaly constant coefficient (float)
    :param a1: Anomaly multiplicative coefficient (float)
    :param a1multiplier: Either 'salinity' to use salinity as the multiplier,
        or 'sigma-theta' to use density (sigma-theta) as the
        multiplier. Defaults to 'salinity'.

    :return: Potential temperature anomaly in ITS-90 degrees C (ndarray)
    """

    # Calculate anomaly correction if coefficients are non-zero
    if a0 != 0.0 or a1 != 0.0:
        po_temp_90_c = sw.ptmp(salinity, temperature, pressure, 0)
        if a1multiplier == "sigma-theta":
            # TODO: should we
            density_ref = sw.pden(salinity, po_temp_90_c, pressure, 0)
            anomaly = po_temp_90_c - (a0 + a1 * density_ref)
        else:
            # Use practical salinity (default)
            anomaly = po_temp_90_c - (a0 + a1 * salinity)
    else:
        anomaly = np.full_like(temperature, const.FLAG_VALUE)  # No anomaly correction needed

    # Potential temperature anomaly is the anomaly correction applied
    return anomaly


def derive_thermosteric_anomaly(
    salinity: np.ndarray,
    temperature: np.ndarray,
) -> np.ndarray:
    """Derive thermosteric anomaly from salinity and temperature using EOS-80.

    Calculates the specific volume anomaly using EOS-80 formulations.

    :param salinity: Measured salinity in practical salinity units (PSU)
    :param temperature: Temperature in ITS-90 degrees C

    :return: Thermosteric anomaly values
    """

    density_val = sw.dens0(salinity, temperature) - 1000
    return 1.0e5 * ((1000.0 / (1000.0 + density_val)) - 0.97266)


def derive_sound_velocity_c(
    salinity: np.ndarray,
    temperature: np.ndarray,
    pressure: np.ndarray,
) -> np.ndarray:
    """Derive sound velocity using Chen and Millero (1977) formula with EOS-80.

    Calculates sound velocity in seawater from salinity, temperature,
    and pressure using the Chen and Millero algorithm based on EOS-80.

    :param salinity: Measured salinity in practical salinity units (PSU)
    :param temperature: Temperature in ITS-90 degrees C
    :param pressure: Measured pressure in decibars

    :return: Sound velocity in m/s
    """

    s_arr, t_arr, p_arr = np.broadcast_arrays(salinity, temperature, pressure)

    s = s_arr.astype(float, copy=False)
    t = t_arr.astype(float, copy=False)
    p = p_arr.astype(float, copy=False) / 10.0

    s = np.maximum(s, 0.0)
    sr = np.sqrt(s)

    d = 1.727e-3 - 7.9836e-6 * p
    b1 = 7.3637e-5 + 1.7945e-7 * t
    b0 = -1.922e-2 - 4.42e-5 * t
    b = b0 + b1 * p

    a3 = (-3.389e-13 * t + 6.649e-12) * t + 1.100e-10
    a2 = ((7.988e-12 * t - 1.6002e-10) * t + 9.1041e-9) * t - 3.9064e-7
    a1 = (((-2.0122e-10 * t + 1.0507e-8) * t - 6.4885e-8) * t - 1.2580e-5) * t + 9.4742e-5
    a0 = (((-3.21e-8 * t + 2.006e-6) * t + 7.164e-5) * t - 1.262e-2) * t + 1.389
    a = ((a3 * p + a2) * p + a1) * p + a0

    c3 = (-2.3643e-12 * t + 3.8504e-10) * t - 9.7729e-9
    c2 = (((1.0405e-12 * t - 2.5335e-10) * t + 2.5974e-8) * t - 1.7107e-6) * t + 3.1260e-5
    c1 = (((-6.1185e-10 * t + 1.3621e-7) * t - 8.1788e-6) * t + 6.8982e-4) * t + 0.153563
    c0 = (
        (((3.1464e-9 * t - 1.47800e-6) * t + 3.3420e-4) * t - 5.80852e-2) * t + 5.03711
    ) * t + 1402.388
    c = ((c3 * p + c2) * p + c1) * p + c0

    return c + (a + b * sr + d * s) * s


def derive_sound_velocity_d(
    salinity: np.ndarray,
    temperature: np.ndarray,
    pressure: np.ndarray,
) -> np.ndarray:
    """Derive sound velocity using Del Grosso (1974) formula with EOS-80.

    Calculates sound velocity in seawater from salinity, temperature,
    and pressure using the Del Grosso algorithm based on EOS-80.

    :param salinity: Measured salinity in practical salinity units (PSU)
    :param temperature: Temperature in ITS-90 degrees C
    :param pressure: Measured pressure in decibars

    :return: Sound velocity in m/s
    """

    s_arr, t_arr, p_arr = np.broadcast_arrays(salinity, temperature, pressure)

    s = s_arr.astype(float, copy=False)
    t = t_arr.astype(float, copy=False)
    p = p_arr.astype(float, copy=False) / 9.80665

    c000 = 1402.392

    dct = (0.501109398873e1 - (0.550946843172e-1 - 0.22153596924e-3 * t) * t) * t
    dcs = (0.132952290781e1 + 0.128955756844e-3 * s) * s
    dcp = (0.156059257041e0 + (0.244998688441e-4 - 0.83392332513e-8 * p) * p) * p
    dcstp = (
        -0.127562783426e-1 * t * s
        + 0.635191613389e-2 * t * p
        + 0.265484716608e-7 * t * t * p * p
        - 0.159349479045e-5 * t * p * p
        + 0.522116437235e-9 * t * p * p * p
        - 0.438031096213e-6 * t * t * t * p
        - 0.161674495909e-8 * s * s * p * p
        + 0.968403156410e-4 * t * t * s
        + 0.485639620015e-5 * t * s * s * p
        - 0.340597039004e-3 * t * s * p
    )

    return c000 + dct + dcs + dcp + dcstp


def derive_sound_velocity_w(
    salinity: np.ndarray,
    temperature: np.ndarray,
    pressure: np.ndarray,
) -> np.ndarray:
    """Derive sound velocity using Wilson (1960) formula with EOS-80.

    Calculates sound velocity in seawater from salinity, temperature,
    and pressure using the Wilson algorithm based on EOS-80.

    :param salinity: Measured salinity in practical salinity units (PSU)
    :param temperature: Temperature in ITS-90 degrees C
    :param pressure: Measured pressure in decibars

    :return: Sound velocity in m/s
    """

    s_arr, t_arr, p_arr = np.broadcast_arrays(salinity, temperature, pressure)

    s = s_arr.astype(float, copy=False)
    t = t_arr.astype(float, copy=False)
    p = p_arr.astype(float, copy=False)

    pr = 0.1019716 * (p + 10.1325)
    sd = s - 35.0

    a = (((7.9851e-6 * t - 2.6045e-4) * t - 4.4532e-2) * t + 4.5721) * t + 1449.14
    sv = (7.7711e-7 * t - 1.1244e-2) * t + 1.39799
    v0 = (1.69202e-3 * sd + sv) * sd + a

    a = ((4.5283e-8 * t + 7.4812e-6) * t - 1.8607e-4) * t + 0.16072
    sv = (1.579e-9 * t + 3.158e-8) * t + 7.7016e-5
    v1 = sv * sd + a

    a = (1.8563e-9 * t - 2.5294e-7) * t + 1.0268e-5
    sv = -1.2943e-7 * sd + a
    a = -1.9646e-10 * t + 3.5216e-9

    return (((-3.3603e-12 * pr + a) * pr + sv) * pr + v1) * pr + v0


def derive_specific_conductance(
    temperature: np.ndarray,
    conductivity: np.ndarray,
    to_units: Literal["uS/cm", "umhos/cm", "mS/cm", "mmhos/cm"],
) -> np.ndarray:
    """Derive specific conductance from temperature and conductivity using EOS-80.

    Calculates specific conductance using temperature and conductivity values
    with EOS-80 formulations.
    Note: 'uS/cm' and 'umhos/cm' are equivalent, as are 'mS/cm' and 'mmhos/cm'.
    All equation information comes from SBE Data Processing code.
    It is intentional that 'uS/cm' and 'umhos/cm' map to the same value, same with 'mS/cm' and 'mmhos/cm'

    :param temperature: Temperature values in ITS-90 degrees C
    :param conductivity: Conductivity values in S/m
    :param to_units: Target units for output (uS/cm, umhos/cm, mS/cm, or mmhos/cm)

    :return: Specific conductance in the specified units
    """
    divide_by = (temperature - 25.0) * const.THERMAL_CONDUCTIVITY_COEFF + 1.0
    if to_units == "uS/cm" or to_units == "umhos/cm":
        return 10000.0 * conductivity / divide_by
    elif to_units == "mmhos/cm" or to_units == "mS/cm":
        return 10.0 * conductivity / divide_by
    else:
        return np.full(len(conductivity), const.FLAG_VALUE)


def derive_average_sound_velocity(
    depth_m: np.ndarray,
    pressure: np.ndarray,
    salinity: np.ndarray,
    sound_velocity_mps: np.ndarray,
    min_pressure: float = -np.inf,
    min_salinity: float = -np.inf,
    delta_d_min: float = 0.5,
) -> np.ndarray:
    """Compute Average Sound Velocity (ASV) over a profile using EOS-80.

    Data is expected to be binned by depth prior to running this function,
    but utilizes delta_d_min (hardcoded to 0.5 m in SBE Data Processing) to determine when to update the ASV value.
    Uses EOS-80 formulations for sound velocity calculations.

    :param depth_m: Depth in meters
    :param pressure: Pressure series (used for validity and min_pressure threshold)
    :param salinity: Salinity series (used for min_salinity threshold)
    :param sound_velocity_mps: Sound velocity in m/s (precomputed)
    :param min_pressure: Ignore rows where pressure < min_pressure, defaults to -inf
    :param min_salinity: Ignore rows where salinity < min_salinity, defaults to -inf
    :param delta_d_min: Minimum depth increment (m) from the last accepted depth to update ASV, defaults to 0.5

    :return: ASV in m/s at each sample. NaN before first valid ASV point
    """
    depth_m = np.asarray(depth_m, dtype=np.float64)
    pressure = np.asarray(pressure, dtype=np.float64)
    salinity = np.asarray(salinity, dtype=np.float64)
    sound_velocity_mps = np.asarray(sound_velocity_mps, dtype=np.float64)

    if not (depth_m.ndim == pressure.ndim == salinity.ndim == sound_velocity_mps.ndim == 1):
        raise ValueError("All inputs must be 1D arrays.")
    n = depth_m.size
    if not (pressure.size == salinity.size == sound_velocity_mps.size == n):
        raise ValueError("All inputs must have the same length.")

    out = np.full(n, np.nan, dtype=np.float64)

    valid = (sound_velocity_mps != 0.0) & (pressure >= min_pressure) & (salinity >= min_salinity)

    valid_idx = np.flatnonzero(valid)
    if valid_idx.size == 0:
        return out

    first = int(valid_idx[0])
    d0 = depth_m[first]
    sv0 = sound_velocity_mps[first]
    dsum = d0
    tsum = d0 / sv0
    asv = sv0

    out[first] = asv
    if first > 0:
        out[:first] = np.nan

    # Hold the current ASV value on each subsequent row, updating only when
    # the depth increment threshold is met.
    for i in range(first + 1, n):
        if valid[i]:
            delta_d = depth_m[i] - d0
            if delta_d >= delta_d_min:
                tsum += delta_d / sound_velocity_mps[i]
                dsum += delta_d
                asv = dsum / tsum
                d0 = depth_m[i]
        out[i] = asv

    return out


def derive_gpa(
    sva: np.ndarray,
    pressure: np.ndarray,
) -> np.ndarray:
    """Derive Geopotential Anomaly from specific volume anomaly and pressure using EOS-80.

    This is a vectorized port of the legacy C++ ``CompGPA`` logic:
    - For each scan, compute the change in GPA based on SVA and pressure changes
    - GPA is accumulated when pressure increases
    - Uses the average SVA between consecutive points and the pressure delta
    - Based on EOS-80 specific volume anomaly calculations

    :param sva: specific volume anomaly values
    :param pressure: pressure values in dbar

    :return: geopotential anomaly values
    """
    sva_arr, pressure_arr = np.broadcast_arrays(sva, pressure)

    sva = sva_arr.astype(float, copy=False)
    p = pressure_arr.astype(float, copy=False)

    # Compute pressure differences (0 for first element)
    p_diff = np.diff(p, prepend=p[0])

    # Get previous SVA values (0 for first element)
    sva_prev = np.concatenate(([0.0], sva[:-1]))

    # Average SVA between consecutive points
    sva_avg = (sva_prev + sva) / 2.0

    # Compute delta GPA (only when pressure increases)
    delta_gpa = np.where(p_diff > 0, sva_avg * (p_diff / 10000.0), 0.0)

    # Cumulative sum to get GPA
    gpa = np.cumsum(delta_gpa)

    return gpa
