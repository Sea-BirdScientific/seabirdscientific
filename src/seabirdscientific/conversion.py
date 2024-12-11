#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""A collection of raw data conversion functions.

Functions:

    convert_temperature (np.ndarray, TemperatureCoefficients, str, str, bool)
    convert_pressure ( np.ndarray, np.ndarray, PressureCoefficients, str)
    convert_conductivity (np.ndarray, np.ndarray, np.ndarray, ConductivityCoefficients)
    potential_density_from_t_s_p (np.ndarray, np.ndarray, np.ndarray, float, float, float)
    potential_density_from_t_c_p (np.ndarray, np.ndarray, np.ndarray, float, float, float)
    density_from_t_s_p (np.ndarray, np.ndarray, np.ndarray, float, float)
    density_from_t_c_p (np.ndarray, np.ndarray, np.ndarray, float, float)
    depth_from_pressure (np.ndarray, float, depth_units="m", pressure_units="dbar")
    convert_sbe63_oxygen (np.ndarray, np.ndarray, np.ndarray, np.ndarray, Oxygen63Coefficients, Thermistor63Coefficients, str)
    convert_sbe63_thermistor (np.ndarray, Thermistor63Coefficients)
    convert_sbe43_oxygen (np.ndarray, np.ndarray, np.ndarray, np.ndarray, Oxygen43Coefficients, bool, bool, float, float)
    convert_oxygen_to_mg_per_l (np.ndarray)
    convert_oxygen_to_umol_per_kg (np.ndarray, np.ndarray)
    convert_eco_chlorophylla_val (float, ChlorophyllACoefficients)
    convert_eco_turbidity_val (float, TurbidityCoefficients)
    convert_sbe18_ph_val(float, float, PH18Coefficients)
    convert_par_logarithmic_val(float, PARCoefficients)

"""

# Native imports
from math import sqrt, e, log, exp, floor
from typing import Literal

# Third-party imports
import gsw
import numpy as np
from scipy import stats

# Sea-Bird imports

# Internal imports
from .cal_coefficients import (
    ChlorophyllACoefficients,
    ConductivityCoefficients,
    Oxygen43Coefficients,
    Oxygen63Coefficients,
    PARCoefficients,
    PH18Coefficients,
    PressureCoefficients,
    TemperatureCoefficients,
    Thermistor63Coefficients,
    TurbidityCoefficients,
)


DBAR_TO_PSI = 1.450377
PSI_TO_DBAR = 0.6894759
OXYGEN_PHASE_TO_VOLTS = 39.457071
KELVIN_OFFSET_0C = 273.15
KELVIN_OFFSET_25C = 298.15
OXYGEN_MLPERL_TO_MGPERL = 1.42903
OXYGEN_MLPERL_TO_UMOLPERKG = 44660
ITS90_TO_IPTS68 = 1.00024  # taken from https://blog.seabird.com/ufaqs/what-is-the-difference-in-temperature-expressions-between-ipts-68-and-its-90/


def convert_temperature(
    temperature_counts_in: np.ndarray,
    coefs: TemperatureCoefficients,
    standard: Literal["ITS90", "IPTS68"] = "ITS90",
    units: Literal["C", "F"] = "C",
    use_mv_r: bool = False,
):
    """Returns the value after converting it to degrees C, ITS-90.

    Data is expected to be raw data from an instrument in A/D counts

    Args:
        temperature_counts_in (np.ndarray): temperature value to convert in A/D counts
        coefs (TemperatureCoefficients) calibration coefficients for the temperature sensor
        standard (str): whether to use ITS90 or to use IPTS-68 calibration standard
        units (str): whether to use celsius or to convert to fahrenheit
        use_mv_r (bool): true to perform extra conversion steps required
        by some instruments (check the cal sheet to see if this is required)

    Returns:
        int: temperature val converted to ITS-90 degrees C
    """

    if use_mv_r:
        mv = (temperature_counts_in - 524288) / 1.6e007
        r = (mv * 2.900e009 + 1.024e008) / (2.048e004 - mv * 2.0e005)
        temperature_counts = r
    else:
        temperature_counts = temperature_counts_in

    log_t = np.log(temperature_counts)
    temperature = (
        1 / (coefs.a0 + coefs.a1 * log_t + coefs.a2 * log_t**2 + coefs.a3 * log_t**3)
    ) - KELVIN_OFFSET_0C

    if standard == "IPTS68":
        temperature *= ITS90_TO_IPTS68
    if units == "F":
        temperature = temperature * 9 / 5 + 32  # Convert C to F

    return temperature


def convert_pressure(
    pressure_count: np.ndarray,
    compensation_voltage: np.ndarray,
    coefs: PressureCoefficients,
    units: Literal["dbar", "psia"] = "psia",
):
    """Returns the value after converting it to psia (pounds per square inch, abolute).

    pressure_count and compensation_voltage are expected to be raw data from an instrument in A/D counts

    Args:
        pressure_count (np.ndarray): pressure value to convert, in A/D counts
        compensation_voltage (np.ndarray): pressure temperature compensation voltage,
            in counts or volts depending on the instrument
        coefs (PressureCoefficients): calibration coefficients for the pressure sensor
        units (str): whether or not to use psia or dbar as the returned unit type

    Returns:
        int: pressure val in PSIA
    """
    sea_level_pressure = 14.7

    t = (
        coefs.ptempa0
        + coefs.ptempa1 * compensation_voltage
        + coefs.ptempa2 * compensation_voltage**2
    )
    x = pressure_count - coefs.ptca0 - coefs.ptca1 * t - coefs.ptca2 * t**2
    n = x * coefs.ptcb0 / (coefs.ptcb0 + coefs.ptcb1 * t + coefs.ptcb2 * t**2)
    pressure = coefs.pa0 + coefs.pa1 * n + coefs.pa2 * n**2 - sea_level_pressure

    if units == "dbar":
        pressure *= PSI_TO_DBAR

    return pressure


def convert_conductivity(
    conductivity_count: np.ndarray,
    temperature: np.ndarray,
    pressure: np.ndarray,
    coefs: ConductivityCoefficients,
):
    """Returns the value after converting it to S/m.

    Data is expected to be raw data from instrument in A/D counts

    Args:
        conductivity_count (np.ndarray): conductivity value to convert, in A/D counts
        temperature (np.ndarray): temperature value to use are reference, in degrees C
        pressure (np.ndarray): pressure value to use are reference, in dbar
        coefs (float): coefs calibration coefficient for the conductivity sensor
    Returns:
        Decimal: conductivity val converted to S/m
    """
    f = conductivity_count * np.sqrt(1 + coefs.wbotc * temperature) / 1000
    numerator = coefs.g + coefs.h * f**2 + coefs.i * f**3 + coefs.j * f**4
    denominator = 1 + coefs.ctcor * temperature + coefs.cpcor * pressure
    return numerator / denominator


def potential_density_from_t_s_p(
    temperature: np.ndarray,
    salinity: np.ndarray,
    pressure: np.ndarray,
    lon=0.0,
    lat=0.0,
    reference_pressure=0.0,
):
    """Derive potential density from measured temperature, salinity, and pressure.

    See: TEOS_10.cpp line 953

    Args:
        temperature (np.ndarray): Measure temperature in degrees C
        salinity (np.ndarray): Measured salinity in practical salinity units
        pressure (np.ndarray): Measured pressure in decibars
        lon (float): Longitude
        lat (float): Latitude
        reference_pressure (float, optional): Reference pressure in decibars. Defaults to 0.0.

    Returns:
        np.ndarray: Potential density in kg/m^3
    """

    absolute_salinity = gsw.SA_from_SP(salinity, pressure, lon, lat)
    conservative_temperature = gsw.CT_from_t(absolute_salinity, temperature, pressure)
    potential_density = (
        gsw.rho(absolute_salinity, conservative_temperature, reference_pressure) - 1000
    )
    return potential_density


def potential_density_from_t_c_p(
    temperature: np.ndarray,
    conductivity: np.ndarray,
    pressure: np.ndarray,
    lon=0.0,
    lat=0.0,
    reference_pressure=0.0,
):
    """Derive potential density from measured temperature, salinity, and pressure.

    See: TEOS_10.cpp line 953

    Args:
        temperature (np.ndarray): Measure temperature in degrees C
        conductivity (np.ndarray): Measured conductivity in mSiemens/cm
        pressure (np.ndarray): Measured pressure in decibars
        lon (float): Longitude
        lat (float): Latitude
        reference_pressure (float, optional): Reference pressure in decibars. Defaults to 0.0.

    Returns:
        np.ndarray: Potential density in kg/m^3
    """

    salinity = gsw.SP_from_C(conductivity, temperature, pressure)
    return potential_density_from_t_s_p(
        temperature, salinity, pressure, lon, lat, reference_pressure
    )


def density_from_t_s_p(
    temperature: np.ndarray,
    salinity: np.ndarray,
    pressure: np.ndarray,
    lon=0.0,
    lat=0.0,
):
    """Derive potential density from measured temperature, salinity, and pressure.

    See: TEOS_10.cpp line 953

    Args:
        temperature (np.ndarray): Measure temperature in degrees C
        salinity (np.ndarray): Measured salinity in practical salinity units
        pressure (np.ndarray): Measured pressure in decibars
        lon (float): Longitude
        lat (float): Latitude

    Returns:
        np.ndarray: Potential density in kg/m^3
    """

    absolute_salinity = gsw.SA_from_SP(salinity, pressure, lon, lat)
    conservative_temperature = gsw.CT_from_t(absolute_salinity, temperature, pressure)
    density = gsw.rho(absolute_salinity, conservative_temperature, pressure)
    return density


def density_from_t_c_p(
    temperature: np.ndarray,
    conductivity: np.ndarray,
    pressure: np.ndarray,
    lon=0.0,
    lat=0.0,
):
    """Derive potential density from measured temperature, salinity, and pressure.

    See: TEOS_10.cpp line 953

    Args:
        temperature (np.ndarray): Measure temperature in degrees C
        conductivity (np.ndarray): Measured conductivity in mSiemens/cm
        pressure (np.ndarray): Measured pressure in decibars
        lon (float): Longitude
        lat (float): Latitude

    Returns:
        np.ndarray: Potential density in kg/m^3
    """

    salinity = gsw.SP_from_C(conductivity, temperature, pressure)
    return density_from_t_s_p(temperature, salinity, pressure, lon, lat)


def depth_from_pressure(
    pressure_in: np.ndarray,
    latitude: float,
    depth_units: Literal["m", "ft"] = "m",
    pressure_units: Literal["dbar", "psi"] = "dbar",
):
    """Derive depth from pressure and latitude.

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


def convert_sbe63_oxygen(
    raw_oxygen_phase: np.ndarray,
    thermistor: np.ndarray,
    pressure: np.ndarray,
    salinity: np.ndarray,
    coefs: Oxygen63Coefficients,
    thermistor_coefs: Thermistor63Coefficients,
    thermistor_units: Literal["volts", "C"] = "volts",
):
    """Returns the data after converting it to ml/l.

        raw_oxygen_phase is expected to be in raw phase, raw_thermistor_temp in counts, pressure in dbar, and salinity in practical salinity (PSU)

        Args:
            raw_oxygen_phase (np.ndarray): SBE63 phase value, in microseconds
    #         raw_thermistor_temp (np.ndarray): SBE63 thermistor data to use are reference, in counts
            pressure (np.ndarray): Converted pressure value from the attached CTD, in dbar
            salinity (np.ndarray): Converted salinity value from the attached CTD, in practical salinity PSU
            coefs (Oxygen63Coefficients): calibration coefficients for the SBE63 sensor
        Returns:
            np.ndarray: converted Oxygen value, in ml/l
    """
    if thermistor_units == "volts":
        temperature = convert_sbe63_thermistor(thermistor, thermistor_coefs)
    elif thermistor_units == "C":
        temperature = thermistor

    oxygen_volts = raw_oxygen_phase / OXYGEN_PHASE_TO_VOLTS  # from the manual
    # O2 (ml/L) = [((a0 + a1T + a2(V^2)) / (b0 + b1V) – 1) / Ksv] [SCorr] [PCorr]

    # Ksv = c0 + c1T + c2 (T^2)
    ksv = coefs.c0 + coefs.c1 * temperature + coefs.c2 * temperature**2

    # SCorr = exp [S * (SolB0 + SolB1 * Ts + SolB2 * Ts^2 + SolB3 * Ts^3) + SolC0 * S^2]
    # The following correction coefficients are all constants
    SOL_B0 = -6.24523e-3
    SOL_B1 = -7.37614e-3
    SOL_B2 = -1.0341e-2
    SOL_B3 = -8.17083e-3
    SOL_C0 = -4.88682e-7

    # Ts = ln [(298.15 – T) / (273.15 + T)]
    ts = np.log((KELVIN_OFFSET_25C - temperature) / (KELVIN_OFFSET_0C + temperature))
    s_corr_exp = (
        salinity * (SOL_B0 + SOL_B1 * ts + SOL_B2 * ts**2 + SOL_B3 * ts**3) + SOL_C0 * salinity**2
    )
    s_corr = e**s_corr_exp

    # Pcorr = exp (E * P / K)
    K = temperature + KELVIN_OFFSET_0C
    # temperature in Kelvin
    p_corr_exp = (e * pressure) / K
    p_corr = e**p_corr_exp

    ox_val = (
        (
            (
                (coefs.a0 + coefs.a1 * temperature + coefs.a2 * oxygen_volts**2)
                / (coefs.b0 + coefs.b1 * oxygen_volts)
                - 1.0
            )
            / ksv
        )
        * s_corr
        * p_corr
    )

    return ox_val


def convert_sbe63_thermistor(
    instrument_output: np.ndarray,
    coefs: Thermistor63Coefficients,
):
    """Converts a SBE63 thermistor raw output array to temperature in ITS-90 deg C.

    Args:
        instrument_output (np.ndarray) raw values from the thermistor
        coefs (Thermisto63Coefficients): calibration coefficients for the thermistor in the SBE63 sensor

    Returns:
        np.ndarray: converted thermistor temperature values in ITS-90 deg C
    """
    logVal = np.log((100000 * instrument_output) / (3.3 - instrument_output))
    temperature = (
        1 / (coefs.ta0 + coefs.ta1 * logVal + coefs.ta2 * logVal**2 + coefs.ta3 * logVal**3)
        - KELVIN_OFFSET_0C
    )
    return temperature


def convert_sbe43_oxygen(
    voltage: np.ndarray,
    temperature: np.ndarray,
    pressure: np.ndarray,
    salinity: np.ndarray,
    coefs: Oxygen43Coefficients,
    apply_tau_correction: bool = False,
    apply_hysteresis_correction: bool = False,
    window_size: float = 1,
    sample_interval: float = 1,
):
    """Returns the data after converting it to ml/l.

    voltage is expected to be in volts, temperature in deg c, pressure in dbar, and salinity in practical salinity (PSU)
    All equation information comes from the June 2013 revision of the SBE43 manual

    Args:
        voltage (float): SBE43 voltage
        temperature (float): temperature value converted to deg C
        pressure (float): Converted pressure value from the attached CTD, in dbar
        salinity (float): Converted salinity value from the attached CTD, in practical salinity PSU
        coefs (Oxygen43Coefficients): calibration coefficients for the SBE43 sensor
        apply_tau_correction (bool): whether or not to run tau correction
        apply_hysteresis_correction (bool): whether or not to run hysteresis correction
        window_size (float): size of the window to use for tau correction, if applicable. In seconds.
        sample_interval (float): sample rate of the data to be used for tau correction, if applicable. In seconds.

    Returns:
        np.ndarray: converted Oxygen values, in ml/l
    """
    # start with all 0 for the dvdt
    dvdt_values = np.zeros(len(voltage))
    if apply_tau_correction:
        # Calculates how many scans to have on either side of our median point, accounting for going out of index bounds
        scans_per_side = floor(window_size / 2 / sample_interval)
        for i in range(scans_per_side, len(voltage) - scans_per_side):
            ox_subset = voltage[i - scans_per_side : i + scans_per_side + 1]

            time_subset = np.arange(
                0, len(ox_subset) * sample_interval, sample_interval, dtype=float
            )

            result = stats.linregress(time_subset, ox_subset)

            dvdt_values[i] = result.slope

    correct_ox_voltages = voltage.copy()
    if apply_hysteresis_correction:
        # Hysteresis starts at 1 because 0 can't be corrected
        for i in range(1, len(correct_ox_voltages)):
            # All Equation info from APPLICATION NOTE NO. 64-3
            d = 1 + coefs.h1 * (np.exp(pressure[i] / coefs.h2) - 1)
            c = np.exp(-1 * sample_interval / coefs.h3)
            ox_volts = correct_ox_voltages[i] + coefs.v_offset

            prev_ox_volts_new = correct_ox_voltages[i - 1] + coefs.v_offset
            ox_volts_new = ((ox_volts + prev_ox_volts_new * c * d) - (prev_ox_volts_new * c)) / d
            ox_volts_final = ox_volts_new - coefs.v_offset
            correct_ox_voltages[i] = ox_volts_final

    oxygen = _convert_sbe43_oxygen(
        correct_ox_voltages,
        temperature,
        pressure,
        salinity,
        coefs,
        dvdt_values,
    )
    return oxygen


def _convert_sbe43_oxygen(
    voltage: np.ndarray,
    temperature: np.ndarray,
    pressure: np.ndarray,
    salinity: np.ndarray,
    coefs: Oxygen43Coefficients,
    dvdt_value: np.ndarray,
):
    """Returns the data after converting it to ml/l.

    voltage is expected to be in volts, temperature in deg c, pressure in dbar, and salinity in practical salinity (PSU)
    All equation information comes from the June 2013 revision of the SBE43 manual.
    Expects that hysteresis correction is already performed on the incoming voltage, if desired.

    Args:
        voltage (float): SBE43 voltage
        temperature (float): temperature value converted to deg C
        pressure (float): Converted pressure value from the attached CTD, in dbar
        salinity (float): Converted salinity value from the attached CTD, in practical salinity PSU
        coefs (Oxygen43Coefficients): calibration coefficients for the SBE43 sensor
        dvdt_value (float): derivative value of voltage with respect to time at this point. Expected to be 0 if not using Tau correction

    Returns:
        float: converted Oxygen value, in ml/l
    """

    # Oxygen Solubility equation constants, From SBE43 Manual Appendix A
    A0 = 2.00907
    A1 = 3.22014
    A2 = 4.0501
    A3 = 4.94457
    A4 = -0.256847
    A5 = 3.88767
    B0 = -0.00624523
    B1 = -0.00737614
    B2 = -0.010341
    B3 = -0.00817083
    C0 = -0.000000488682

    ts = np.log((KELVIN_OFFSET_25C - temperature) / (KELVIN_OFFSET_0C + temperature))
    a_term = A0 + A1 * ts + A2 * ts**2 + A3 * ts**3 + A4 * ts**4 + A5 * ts**5
    b_term = salinity * (B0 + B1 * ts + B2 * ts**2 + B3 * ts**3)
    c_term = C0 * salinity**2
    solubility = np.exp(a_term + b_term + c_term)

    # Tau correction
    tau = coefs.tau_20 * np.exp(coefs.d1 * pressure + coefs.d2 * (temperature - 20)) * dvdt_value

    soc_term = coefs.soc * (voltage + coefs.v_offset + tau)
    temp_term = 1.0 + coefs.a * temperature + coefs.b * temperature**2 + coefs.c * temperature**3
    oxygen = (
        soc_term
        * solubility
        * temp_term
        * np.exp((coefs.e * pressure) / (temperature + KELVIN_OFFSET_0C))
    )
    return oxygen


def convert_oxygen_to_mg_per_l(ox_values: np.ndarray):
    """Converts given oxygen values to milligrams/Liter.

    Expects oxygen values to be in Ml/L

    Args:
        ox_values (np.ndarray): oxygen values, already converted to ml/L

    Returns:
        np.ndarray: oxygen values converted to milligrams/Liter
    """
    # From the SBE43 and SBE63 manual
    return ox_values * OXYGEN_MLPERL_TO_MGPERL


def convert_oxygen_to_umol_per_kg(ox_values: np.ndarray, potential_density: np.ndarray):
    """Converts given oxygen values to milligrams/Liter.

    Expects oxygen values to be in Ml/L
    Note: Sigma-Theta is expected to be calculated via gsw_sigma0, meaning is it technically potential density anomaly.
    Calculating using gsw_rho(SA, CT, p_ref = 0) results in actual potential density, but this function already does the converison,
    So values will need to have 1000 subtracted from them before being passed into this function.
    The function is done this way to stay matching to the manual for the SBE63 and SBE43, but the results of either method are identical.

    Args:
        ox_values (np.ndarray): oxygen values, already converted to ml/L
        potential_density (np.ndarray): potential density (sigma-theta) values.
                                        Expected to be the same length as ox_values

    Returns:
        np.ndarray: oxygen values converted to milligrams/Liter
    """
    # From the SBE43 and SBE63 manual
    convertedVals = np.divide(ox_values * OXYGEN_MLPERL_TO_UMOLPERKG, (potential_density + 1000))
    return convertedVals


def convert_eco_chlorophylla(
    raw_chlorophyll_a: float,
    coefs: ChlorophyllACoefficients,
):
    """Converts a raw value for chlorophyll-a channel on a ECO-FLNTU or ECO-FL.

    All equation information comes from ECO-FLNTU calibration sheets

    Args:
        raw_chlorophyll_a (float): raw counts for digital, raw volts for analog
        coefs (ChlorophyllACoefficients): calibration coefficients for clorophyll-a

    Returns:
        float: converted chlorophyll-a in μg/l
    """
    chlorophylla = coefs.scalar * (raw_chlorophyll_a - coefs.v_blank)

    return chlorophylla


def convert_eco_turbidity(
    raw_turbidity: float,
    coefs: TurbidityCoefficients,
):
    """Converts a raw value for turbidity channel on a ECO-FLNTU.

    All equation information comes from ECO-FLNTU calibration sheets

    Args:
        rawTurbidity(float): raw counts for digital, raw volts for analog
        coefs (TurbidityCoefficients): calbration coefficients for turbidity

    Returns:
        float: converted turbidity in nephelometric turbidity units (NTU)
    """
    turbidity = coefs.scalar * (raw_turbidity - coefs.dark_voltage)

    return turbidity


def convert_sbe18_ph(
    raw_ph: float,
    temperature: float,
    coefs: PH18Coefficients,
):
    """Converts a raw voltage value for pH
        All equation information comes from application note 18-1
    Args:
        raw_ph (float): raw output voltage from pH sensor (0-5V)
        temperature (float): temperature value to use for temperature compensation in degrees C
        coefs (PH18Coefficients): slope and offset for the pH sensor
    Returns:
        float: converted pH
    """
    pH = 7 + (raw_ph - coefs.offset) / (
        1.98416e-4 * (temperature + KELVIN_OFFSET_0C) * coefs.slope
    )
    return pH


def convert_par_logarithmic(
    raw_par: float,
    coefs: PARCoefficients,
):
    """Converts a raw voltage value for PAR to µmol photons/m2*s
        All equation information comes from application note 96
    Args:
        rawpH (float): raw output voltage from PAR sensor
        coefs (PARCoefficients): calibration coefficients for the PAR sensor
    Returns:
        float: converted PAR in µmol photons/m2*s
    """
    PAR = coefs.multiplier * coefs.im * 10 ** ((raw_par - coefs.a0) / coefs.a1)

    return PAR
