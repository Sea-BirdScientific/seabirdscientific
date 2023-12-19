#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""TODO: conversion docstring"""

# Native imports
from math import sqrt, e, log

# Third-party imports
import gsw
import numpy as np

# Sea-Bird imports

# Internal imports


DBAR_TO_PSI = 1.450377
PSI_TO_DBAR = 0.6894759
OXYGEN_PHASE_TO_VOLTS = 39.457071
KELVIN_OFFSET = 273.15


def convert_temperature_array(
    temperature_counts: np.ndarray,
    a0: float,
    a1: float,
    a2: float,
    a3: float,
    ITS90: bool,
    celsius: bool,
    use_MV_R: bool,
):
    """Returns the data after converting it to degrees C, ITS-90.
        Data is expected to be raw data from instrument in A/D counts
    Args:
        temperature_counts (np.ndarray): temperature data to convert, in A/D counts
        a0 (float): a0 calibration coefficient for the temperature sensor
        a1 (float): a1 calibration coefficient for the temperature sensor
        a2 (float): a2 calibration coefficient for the temperature sensor
        a3 (float): a3 calibration coefficient for the temperature sensor
        ITS90 (bool): whether to use ITS90 or to use IPTS-68 conventions
        celsius (bool): whether to use celsius or to convert to fahrenheit
        use_MV_R (bool): passed into convert_temperature_val_ITS90_c to perform
            extra conversion steps required by some instruments
    Returns:
        ndarray: temperature values converted to ITS-90 degrees C, in the same order as input
    """
    ipts68_converison = 1.00024  # taken from https://blog.seabird.com/ufaqs/what-is-the-difference-in-temperature-expressions-between-ipts-68-and-its-90/
    convert_vectorized = np.vectorize(
        convert_temperature_val_ITS90_c, excluded=["a0", "a1", "a2", "a3", "use_MV_R"]
    )
    result = convert_vectorized(temperature_counts, a0, a1, a2, a3, use_MV_R)
    if not ITS90:
        result = result * ipts68_converison
    if not celsius:
        result = result * 9 / 5 + 32  # Convert C to F
    return result


def convert_temperature_val_ITS90_c(
    temperature_counts_in: int,
    a0: float,
    a1: float,
    a2: float,
    a3: float,
    use_MV_R: bool,
):
    """Returns the value after converting it to degrees C, ITS-90.
        Data is expected to be raw data from instrument in A/D counts
    Args:
        temperature_counts_in (int): temperature value to convert, in A/D counts
        a0 (float): a0 calibration coefficient for the temperature sensor
        a1 (float): a1 calibration coefficient for the temperature sensor
        a2 (float): a2 calibration coefficient for the temperature sensor
        a3 (float): a3 calibration coefficient for the temperature sensor
        use_MV_R (bool): whether to perform extra conversion steps required by some instruments
    Returns:
        int: temperature val converted to ITS-90 degrees C
    """
    if use_MV_R:
        MV = (temperature_counts_in - 524288) / 1.6e007
        R = (MV * 2.900e009 + 1.024e008) / (2.048e004 - MV * 2.0e005)
        temperature_counts = R
    else:
        temperature_counts = temperature_counts_in

    temperature = (
        1
        / (
            a0
            + a1 * np.log(temperature_counts)
            + a2 * np.log(temperature_counts) ** 2
            + a3 * np.log(temperature_counts) ** 3
        )
    ) - 273.15
    return temperature


def convert_pressure_array(
    pressure_counts: np.ndarray,
    compensation_voltages: np.ndarray,
    is_dbar: bool,
    PA0: float,
    PA1: float,
    PA2: float,
    PTEMPA0: float,
    PTEMPA1: float,
    PTEMPA2: float,
    PTCA0: float,
    PTCA1: float,
    PTCA2: float,
    PTCB0: float,
    PTCB1: float,
    PTCB2: float,
):
    """Calls convert_pressure_val_strain on an array of raw pressure data.
        Data is expected to be raw data from instrument in A/D counts
    Args:
        pressure_counts (np.ndarray): pressure data to convert, in A/D counts
        compensation_voltages (np.ndarray): pressure temperature compensation voltages,
            in counts or volts depending on instrument
        is_dbar (bool): whether or not to use psia as the returned unit type. If false, uses dbar
        PA0 (float): PA0 calibration coefficient for the pressure sensor
        PA1 (float): PA1 calibration coefficient for the pressure sensor
        PA2 (float): PA2 calibration coefficient for the pressure sensor
        PTEMPA0 (float): PTEMPA0 calibration coefficient for the pressure sensor
        PTEMPA1 (float): PTEMPA1 calibration coefficient for the pressure sensor
        PTEMPA2 (float): PTEMPA2 calibration coefficient for the pressure sensor
        PTCA0 (float): PTCA0 calibration coefficient for the pressure sensor
        PTCA1 (float): PTCA1 calibration coefficient for the pressure sensor
        PTCA2 (float): PTCA2 calibration coefficient for the pressure sensor
        PTCB0 (float): PTCB0 calibration coefficient for the pressure sensor
        PTCB1 (float): PTCB1 calibration coefficient for the pressure sensor
        PTCB2 (float): PTCB2 calibration coefficient for the pressure sensor
    Returns:
        ndarray: pressure values
    """

    pressure = np.empty(shape=(pressure_counts.size))
    for i in range(0, pressure_counts.size):
        pressure[i] = convert_pressure_val_strain(
            pressure_counts[i],
            compensation_voltages[i],
            is_dbar,
            PA0,
            PA1,
            PA2,
            PTEMPA0,
            PTEMPA1,
            PTEMPA2,
            PTCA0,
            PTCA1,
            PTCA2,
            PTCB0,
            PTCB1,
            PTCB2,
        )
    return pressure


def convert_pressure_val_strain(
    pressure_count: float,
    compensation_voltage: float,
    is_dbar: bool,
    PA0: float,
    PA1: float,
    PA2: float,
    PTEMPA0: float,
    PTEMPA1: float,
    PTEMPA2: float,
    PTCA0: float,
    PTCA1: float,
    PTCA2: float,
    PTCB0: float,
    PTCB1: float,
    PTCB2: float,
):
    """Returns the value after converting it to PSIA (pounds per square inch, abolute)
        pressure_count and compensation_voltage are expected to be raw data from instrument in A/D counts
    Args:
        pressure_count (int): pressure value to convert, in A/D counts
        compensation_voltage (float): pressure temperature compensation voltage,
            in counts or volts depending on the instrument
        PA0 (float): PA0 calibration coefficient for the pressure sensor
        PA1 (float): PA1 calibration coefficient for the pressure sensor
        PA2 (float): PA2 calibration coefficient for the pressure sensor
        PTEMPA0 (float): PTEMPA0 calibration coefficient for the pressure sensor
        PTEMPA1 (float): PTEMPA1 calibration coefficient for the pressure sensor
        PTEMPA2 (float): PTEMPA2 calibration coefficient for the pressure sensor
        PTCA0 (float): PTCA0 calibration coefficient for the pressure sensor
        PTCA1 (float): PTCA1 calibration coefficient for the pressure sensor
        PTCA2 (float): PTCA2 calibration coefficient for the pressure sensor
        PTCB0 (float): PTCB0 calibration coefficient for the pressure sensor
        PTCB1 (float): PTCB1 calibration coefficient for the pressure sensor
        PTCB2 (float): PTCB2 calibration coefficient for the pressure sensor
    Returns:
        int: pressure val in PSIA
    """
    sea_level_pressure = 14.7

    t = PTEMPA0 + PTEMPA1 * compensation_voltage + PTEMPA2 * compensation_voltage**2
    x = pressure_count - PTCA0 - PTCA1 * t - PTCA2 * t**2
    n = x * PTCB0 / (PTCB0 + PTCB1 * t + PTCB2 * t**2)
    pressure = PA0 + PA1 * n + PA2 * n**2 - sea_level_pressure

    if is_dbar:
        pressure *= PSI_TO_DBAR

    return pressure


def convert_conductivity_array(
    conductivity_counts: np.ndarray,
    temperature: np.ndarray,
    pressure: np.ndarray,
    g: float,
    h: float,
    i: float,
    j: float,
    CPcor: float,
    CTcor: float,
    WBOTC: float,
):
    """Returns the data after converting it to Siemens/meter (S/m)
        cond_data is expected to be in raw counts, temp_data in C, and press_data in dbar
    Args:
        conductivity_counts (np.ndarray): conductivity data to convert, in A/D counts
        temperature (np.ndarray): temperature data to use are reference, in degrees C
        pressure (np.ndarray): pressure data to use are reference, in dbar
        g (float): g calibration coefficient for the conductivity sensor
        h (float): h calibration coefficient for the conductivity sensor
        i (float): i calibration coefficient for the conductivity sensor
        j (float): j calibration coefficient for the conductivity sensor
        CPcor (float): CPcor calibration coefficient for the conductivity sensor
        CTcor (float): CTcor calibration coefficient for the conductivity sensor
        WBOTC (float): Wien bridge oscillator temperature coefficient
            https://sndl.ucmerced.edu/files/MHWG/Sensors_and_Loggers/Manuals/37SMmanual34829.pdf
    Returns:
        ndarray: conductivity values converted to S/m, in the same order as input
    """
    conductivity = np.empty(shape=(conductivity_counts.size))
    for index in range(0, conductivity_counts.size):
        conductivity[index] = convert_conductivity_val(
            conductivity_counts[index],
            temperature[index],
            pressure[index],
            g,
            h,
            i,
            j,
            CPcor,
            CTcor,
            WBOTC,
        )

    return conductivity


def convert_conductivity_val(
    conductivity_count: float,
    temperature: float,
    pressure: float,
    g: float,
    h: float,
    i: float,
    j: float,
    CPcor: float,
    CTcor: float,
    WBOTC,
):
    """Returns the value after converting it to S/m
        Data is expected to be raw data from instrument in A/D counts
    Args:
        conductivity_count (np.ndarray): conductivity value to convert, in A/D counts
        temperature (np.ndarray): temperature value to use are reference, in degrees C
        pressure (np.ndarray): pressure value to use are reference, in dbar
        g (float): g calibration coefficient for the conductivity sensor
        h (float): h calibration coefficient for the conductivity sensor
        i (float): i calibration coefficient for the conductivity sensor
        j (float): j calibration coefficient for the conductivity sensor
        CPcor (float): CPcor calibration coefficient for the conductivity sensor
        CTcor (float): CTcor calibration coefficient for the conductivity sensor
        WBOTC (float): Wien bridge oscillator temperature coefficient
            file:///I:/common/calibration/SBE37/calibrationPDFs/C24682.pdf
    Returns:
        Decimal: conductivity val converted to S/m
    """
    f = conductivity_count * sqrt(1 + WBOTC * temperature) / 1000
    numerator = g + h * f**2 + i * f**3 + j * f**4
    denominator = 1 + CTcor * temperature + CPcor * pressure
    return numerator / denominator


def potential_density_from_t_s_p(
    temperature_C: np.ndarray,
    salinity_PSU: np.ndarray,
    pressure_dbar: np.ndarray,
    lon=0.0,
    lat=0.0,
    reference_pressure=0.0,
):
    """Derive potential density from measured temperature, salinity, and pressure
    See: TEOS_10.cpp line 953

    Args:
        temperature_C (np.ndarray): Measure temperature in degrees C
        salinity_PSU (np.ndarray): Measured salinity in practical salinity units
        pressure_dbar (np.ndarray): Measured pressure in decibars
        lon (float): Longitude
        lat (float): Latitude
        reference_pressure (float, optional): Reference pressure in decibars. Defaults to 0.0.

    Returns:
        np.ndarray: Potential density in kg/m^3
    """

    absolute_salinity = gsw.SA_from_SP(salinity_PSU, pressure_dbar, lon, lat)
    conservative_temperature = gsw.CT_from_t(
        absolute_salinity, temperature_C, pressure_dbar
    )
    potential_density = (
        gsw.rho(absolute_salinity, conservative_temperature, reference_pressure) - 1000
    )
    return potential_density


def potential_density_from_t_c_p(
    temperature_C: np.ndarray,
    conductivity_mScm: np.ndarray,
    pressure_dbar: np.ndarray,
    lon=0.0,
    lat=0.0,
    reference_pressure=0.0,
):
    """Derive potential density from measured temperature, salinity, and pressure
    See: TEOS_10.cpp line 953

    Args:
        temperature_C (np.ndarray): Measure temperature in degrees C
        conductivity_mScm (np.ndarray): Measured conductivity in mSiemens/cm
        pressure_dbar (np.ndarray): Measured pressure in decibars
        lon (float): Longitude
        lat (float): Latitude
        reference_pressure (float, optional): Reference pressure in decibars. Defaults to 0.0.

    Returns:
        np.ndarray: Potential density in kg/m^3
    """

    salinity_PSU = gsw.SP_from_C(conductivity_mScm, temperature_C, pressure_dbar)
    return potential_density_from_t_s_p(
        temperature_C, salinity_PSU, pressure_dbar, lon, lat, reference_pressure
    )


def density_from_t_s_p(
    temperature_C: np.ndarray,
    salinity_PSU: np.ndarray,
    pressure_dbar: np.ndarray,
    lon=0.0,
    lat=0.0,
):
    """Derive potential density from measured temperature, salinity, and pressure
    See: TEOS_10.cpp line 953

    Args:
        temperature_C (np.ndarray): Measure temperature in degrees C
        salinity_PSU (np.ndarray): Measured salinity in practical salinity units
        pressure_dbar (np.ndarray): Measured pressure in decibars
        lon (float): Longitude
        lat (float): Latitude

    Returns:
        np.ndarray: Potential density in kg/m^3
    """

    absolute_salinity = gsw.SA_from_SP(salinity_PSU, pressure_dbar, lon, lat)
    conservative_temperature = gsw.CT_from_t(
        absolute_salinity, temperature_C, pressure_dbar
    )
    density = gsw.rho(absolute_salinity, conservative_temperature, pressure_dbar)
    return density


def density_from_t_c_p(
    temperature_C: np.ndarray,
    conductivity_mScm: np.ndarray,
    pressure_dbar: np.ndarray,
    lon=0.0,
    lat=0.0,
):
    """Derive potential density from measured temperature, salinity, and pressure
    See: TEOS_10.cpp line 953

    Args:
        temperature_C (np.ndarray): Measure temperature in degrees C
        conductivity_mScm (np.ndarray): Measured conductivity in mSiemens/cm
        pressure_dbar (np.ndarray): Measured pressure in decibars
        lon (float): Longitude
        lat (float): Latitude

    Returns:
        np.ndarray: Potential density in kg/m^3
    """

    salinity_PSU = gsw.SP_from_C(conductivity_mScm, temperature_C, pressure_dbar)
    return density_from_t_s_p(temperature_C, salinity_PSU, pressure_dbar, lon, lat)


def depth_from_pressure(
    pressure_in: np.ndarray, latitude: float, depth_units="m", pressure_units="dbar"
):
    """Derive depth from pressure and latitude

    Args:
        pressure (np.ndarray): Numpy array of floats representing pressure in dbar or psi
        latitude (float): Latitude (-90.0 to 90.0)
        depth_units (str, optional): 'm' for meters, 'ft' for feet. Defaults to 'm'.
        pressure_units (str, optional): 'dbar' for decibars, 'psi' for PSI. Defaults to 'dbar'.

    Returns:
        np.ndarray: A numpy array representing depth in meters or feet
    """
    pressure = pressure_in.copy()
    if pressure_units == "psi":
        pressure /= DBAR_TO_PSI

    depth = -gsw.z_from_p(pressure, latitude)

    if depth_units == "ft":
        depth *= 3.28084

    return depth


def convert_oxygen_array(
    raw_oxygen_phase: np.ndarray,
    raw_thermistor_temp: np.ndarray,
    pressure: np.ndarray,
    salinity: np.ndarray,
    a0: float,
    a1: float,
    a2: float,
    b0: float,
    b1: float,
    c0: float,
    c1: float,
    c2: float,
    ta0: float,
    ta1: float,
    ta2: float,
    ta3: float,
    e: float,
):
    """Returns the data after converting it to ml/l
        raw_oxygen_phase is expected to be in raw phase, raw_thermistor_temp in counts, pressure in dbar, and salinity in practical salinity (PSU)
    Args:
        raw_oxygen_phase (np.ndarray): SBE63 phase values, in microseconds
        raw_thermistor_temp (np.ndarray): SBE63 thermistor data to use are reference, in counts
        pressure (np.ndarray): Converted pressure values from the attached CTD, in dbar
        salinity (np.ndarraty): Converted salinity values from the attached CTD, in practical salinity PSU
        a0 (float): calibration coefficient for the SBE63 sensor
        a1 (float): calibration coefficient for the SBE63 sensor
        a2 (float): calibration coefficient for the SBE63 sensor
        b0 (float): calibration coefficient for the SBE63 sensor
        b1 (float): calibration coefficient for the SBE63 sensor
        c0 (float): calibration coefficient for the SBE63 sensor
        c1 (float): calibration coefficient for the SBE63 sensor
        c2 (float): calibration coefficient for the SBE63 sensor
        ta0 (float): calibration coefficient for the thermistor in the SBE63 sensor
        ta1 (float): calibration coefficient for the thermistor in the SBE63 sensor
        ta2 (float): calibration coefficient for the thermistor in the SBE63 sensor
        ta3 (float): calibration coefficient for the thermistor in the SBE63 sensor
        e (float): calibration coefficient for the SBE63 sensor
    Returns:
        ndarray: converted Oxygen values, in ml/l, in the same order as input
    """
    thermistor_temperature = convert_SBE63_thermistor_array(
        raw_thermistor_temp, ta0, ta1, ta2, ta3
    )
    # oxygen = np.empty(shape = (raw_oxygen_phase.size))
    convert_vectorized = np.vectorize(
        convert_oxygen_val,
        excluded=["a0", "a1", "a2", "a3", "b0", "b1", "c0", "c1", "c2", "e"],
    )
    oxygen = convert_vectorized(
        raw_oxygen_phase,
        thermistor_temperature,
        pressure,
        salinity,
        a0,
        a1,
        a2,
        b0,
        b1,
        c0,
        c1,
        c2,
        e,
    )

    return oxygen


def convert_oxygen_val(
    raw_oxygen_phase: float,
    temperature: float,
    pressure: float,
    salinity: float,
    a0: float,
    a1: float,
    a2: float,
    b0: float,
    b1: float,
    c0: float,
    c1: float,
    c2: float,
    e: float,
):
    """Returns the data after converting it to ml/l
        raw_oxygen_phase is expected to be in raw phase, raw_thermistor_temp in counts, pressure in dbar, and salinity in practical salinity (PSU)
    Args:
        raw_oxygen_phase (np.ndarray): SBE63 phase value, in microseconds
        temperature (np.ndarray): SBE63 thermistor value converted to deg C
        pressure (np.ndarray): Converted pressure value from the attached CTD, in dbar
        salinity (np.ndarray): Converted salinity value from the attached CTD, in practical salinity PSU
        a0 (float): calibration coefficient for the SBE63 sensor
        a1 (float): calibration coefficient for the SBE63 sensor
        a2 (float): calibration coefficient for the SBE63 sensor
        b0 (float): calibration coefficient for the SBE63 sensor
        b1 (float): calibration coefficient for the SBE63 sensor
        c0 (float): calibration coefficient for the SBE63 sensor
        c1 (float): calibration coefficient for the SBE63 sensor
        c2 (float): calibration coefficient for the SBE63 sensor
        ta0 (float): calibration coefficient for the thermistor in the SBE63 sensor
        ta1 (float): calibration coefficient for the thermistor in the SBE63 sensor
        ta2 (float): calibration coefficient for the thermistor in the SBE63 sensor
        ta3 (float): calibration coefficient for the thermistor in the SBE63 sensor
        e (float): calibration coefficient for the SBE63 sensor
    Returns:
        np.ndarray: converted Oxygen value, in ml/l
    """
    oxygen_volts = raw_oxygen_phase / OXYGEN_PHASE_TO_VOLTS  # from the manual
    # O2 (ml/L) = [((a0 + a1T + a2(V^2)) / (b0 + b1V) – 1) / Ksv] [SCorr] [PCorr]

    # Ksv = c0 + c1T + c2 (T^2)
    ksv = c0 + c1 * temperature + c2 * temperature**2

    # SCorr = exp [S * (SolB0 + SolB1 * Ts + SolB2 * Ts^2 + SolB3 * Ts^3) + SolC0 * S^2]
    # The following correction coefficients are all constants
    Sol_B0 = -6.24523e-3
    Sol_B1 = -7.37614e-3
    Sol_B2 = -1.0341e-2
    Sol_B3 = -8.17083e-3
    Sol_C0 = -4.88682e-7

    # Ts = ln [(298.15 – T) / (273.15 + T)]
    ts = log((298.15 - temperature) / (KELVIN_OFFSET + temperature))
    s_corr_exp = (
        salinity * (Sol_B0 + Sol_B1 * ts + Sol_B2 * ts**2 + Sol_B3 * ts**3)
        + Sol_C0 * salinity**2
    )
    s_corr = e**s_corr_exp

    # Pcorr = exp (E * P / K)
    K = temperature + KELVIN_OFFSET
    # temperature in Kelvin
    p_corr_exp = (e * pressure) / K
    p_corr = e**p_corr_exp

    ox_val = (
        (
            (
                (a0 + a1 * temperature + a2 * oxygen_volts**2)
                / (b0 + b1 * oxygen_volts)
                - 1.0
            )
            / ksv
        )
        * s_corr
        * p_corr
    )

    return ox_val


def convert_SBE63_thermistor_array(
    instrument_output: np.ndarray,
    ta0: float,
    ta1: float,
    ta2: float,
    ta3: float,
):
    """Converts a SBE63 thermistor raw output array to temperature in ITS-90 deg C.
    Args:
        instrument_output (np.ndarray) raw values from the thermistor
        ta0 (float): calibration coefficient for the thermistor in the SBE63 sensor
        ta1 (float): calibration coefficient for the thermistor in the SBE63 sensor
        ta2 (float): calibration coefficient for the thermistor in the SBE63 sensor
        ta3 (float): calibration coefficient for the thermistor in the SBE63 sensor
    Returns:
        np.ndarray: converted thermistor temperature values in ITS-90 deg C
    """
    convert_vectorized = np.vectorize(
        convert_SBE63_thermistor_value, excluded=["ta0", "ta1", "ta2", "ta3"]
    )
    temperature = convert_vectorized(instrument_output, ta0, ta1, ta2, ta3)
    return temperature


def convert_SBE63_thermistor_value(
    instrument_output: float,
    ta0: float,
    ta1: float,
    ta2: float,
    ta3: float,
):
    """Converts a SBE63 thermistor raw output array to temperature in ITS-90 deg C.
        Args:
            instrument_output (np.ndarray) raw values from the thermistor
            ta0 (float): calibration coefficient for the thermistor in the SBE63 sensor
            ta1 (float): calibration coefficient for the thermistor in the SBE63 sensor
            ta2 (float): calibration coefficient for the thermistor in the SBE63 sensor
            ta3 (float): calibration coefficient for the thermistor in the SBE63 sensor
    Returns:
            np.ndarray: converted thermistor temperature values in ITS-90 deg C
    """
    logVal = log((100000 * instrument_output) / (3.3 - instrument_output))
    temperature = (
        1 / (ta0 + ta1 * logVal + ta2 * logVal**2 + ta3 * logVal**3) - KELVIN_OFFSET
    )
    return temperature
