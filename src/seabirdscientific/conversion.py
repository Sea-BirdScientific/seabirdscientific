#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""A collection of raw data conversion functions."""

# Functions:
#     convert_temperature (np.ndarray, TemperatureCoefficients, str, str, bool)
#     convert_pressure ( np.ndarray, np.ndarray, PressureCoefficients, str)
#     convert_conductivity (np.ndarray, np.ndarray, np.ndarray, ConductivityCoefficients)
#     potential_density_from_t_s_p (np.ndarray, np.ndarray, np.ndarray, float, float, float)
#     potential_density_from_t_c_p (np.ndarray, np.ndarray, np.ndarray, float, float, float)
#     density_from_t_s_p (np.ndarray, np.ndarray, np.ndarray, float, float)
#     density_from_t_c_p (np.ndarray, np.ndarray, np.ndarray, float, float)
#     depth_from_pressure (np.ndarray, float, depth_units="m", pressure_units="dbar")
#     convert_sbe63_oxygen (
#         np.ndarray, np.ndarray, np.ndarray, np.ndarray,
#         Oxygen63Coefficients, Thermistor63Coefficients, str
#     )
#     convert_sbe63_thermistor (np.ndarray, Thermistor63Coefficients)
#     convert_sbe43_oxygen (
#         np.ndarray, np.ndarray, np.ndarray, np.ndarray,
#         Oxygen43Coefficients, bool, bool, float, float
#     )
#     convert_oxygen_to_mg_per_l (np.ndarray)
#     convert_oxygen_to_umol_per_kg (np.ndarray, np.ndarray)
#     convert_eco_chlorophylla_val (float, ChlorophyllACoefficients)
#     convert_eco_turbidity_val (float, TurbidityCoefficients)
#     convert_sbe18_ph_val(float, float, PH18Coefficients)
#     convert_par_logarithmic_val(float, PARCoefficients)
#     convert_nitrate(np.ndarray, float, float, str)
#     convert_ph_voltage_counts(np.ndarray)
#     convert_internal_seafet_ph(np.ndarray, temperature: np.ndarray, PHSeaFETInternalCoefficients)
#     convert_external_seafet_ph(
#         np.ndarray, np.ndarray, np.ndarray, np.ndarray,
#         PHSeaFETExternalCoefficients
#     )


# Native imports
from math import e, floor
from typing import Literal

# Third-party imports
import gsw
import numpy as np
from scipy import stats

# Sea-Bird imports

# Internal imports
from .cal_coefficients import (
    ConductivityCoefficients,
    Oxygen43Coefficients,
    Oxygen63Coefficients,
    PARCoefficients,
    PH18Coefficients,
    PHSeaFETInternalCoefficients,
    PHSeaFETExternalCoefficients,
    PressureCoefficients,
    TemperatureCoefficients,
    Thermistor63Coefficients,
    ECOCoefficients,
)


DBAR_TO_PSI = 1.450377
PSI_TO_DBAR = 0.6894759
OXYGEN_PHASE_TO_VOLTS = 39.457071
KELVIN_OFFSET_0C = 273.15
KELVIN_OFFSET_25C = 298.15
OXYGEN_MLPERL_TO_MGPERL = 1.42903
OXYGEN_MLPERL_TO_UMOLPERKG = 44660
# taken from https://blog.seabird.com/ufaqs/what-is-the-difference-in-temperature-expressions-between-ipts-68-and-its-90/ # pylint: disable=line-too-long
ITS90_TO_IPTS68 = 1.00024
# micro moles of nitrate to milligrams of nitrogen per liter
UMNO3_TO_MGNL = 0.014007
# [J K^{-1} mol^{-1}] Gas constant, NIST Reference on Constants retrieved 10-05-2015
R = 8.3144621
# [Coulombs mol^{-1}] Faraday constant, NIST Reference on Constants retrieved 10-05-2015
F = 96485.365


def convert_temperature(
    temperature_counts_in: np.ndarray,
    coefs: TemperatureCoefficients,
    standard: Literal["ITS90", "IPTS68"] = "ITS90",
    units: Literal["C", "F"] = "C",
    use_mv_r: bool = False,
):
    """Returns the value after converting it to degrees C, ITS-90.

    Data is expected to be raw data from an instrument in A/D counts

    :param temperature_counts_in: temperature value to convert in A/D
        counts
    :param coefs: calibration coefficients for the temperature sensor
    :param standard: whether to use ITS90 or to use IPTS-68 calibration
        standard
    :param units: whether to use celsius or to convert to fahrenheit
    :param use_mv_r: true to perform extra conversion steps required by
        some instruments (check the cal sheet to see if this is required)

    :return: temperature val converted to ITS-90 degrees C
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
    """Converts pressure counts to PSIA (pounds per square inch, abolute).

    pressure_count and compensation_voltage are expected to be raw data
    from an instrument in A/D counts

    :param pressure_count: pressure value to convert, in A/D counts
    :param compensation_voltage: pressure temperature compensation
        voltage, in counts or volts depending on the instrument
    :param coefs: calibration coefficients for the pressure sensor
    :param units: whether or not to use psia or dbar as the returned
        unit type

    :return: pressure val in PSIA
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
    """Converts raw conductivity counts to S/m.

    Data is expected to be raw data from instrument in A/D counts

    :param conductivity_count: conductivity value to convert, in A/D
        counts
    :param temperature: reference temperature, in degrees C
    :param pressure: reference pressure, in dbar
    :param coefs: calibration coefficient for the conductivity sensor

    :return: conductivity val converted to S/m
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
    """Derive potential density from measured temperature, salinity, and
    pressure.

    :param temperature: Measure temperature, in degrees C
    :param salinity: Measured salinity, in practical salinity units
    :param pressure: Measured pressure, in decibars
    :param lon: Longitude
    :param lat: Latitude
    :param reference_pressure: Reference pressure in decibars. Defaults
        to 0.0.

    :return: Potential density in kg/m^3
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
    """Derive potential density from measured temperature, salinity, and
    pressure.

    :param temperature: Measure temperature, in degrees C
    :param conductivity: Measured conductivity, in mSiemens/cm
    :param pressure: Measured pressure, in decibars
    :param lon: Longitude
    :param lat: Latitude
    :param reference_pressure: Reference pressure in decibars. Defaults
        to 0.0.

    :return: Potential density in kg/m^3
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
    """Derive potential density from measured temperature, salinity, and
    pressure.

    :param temperature: Measure temperature, in degrees C
    :param salinity: Measured salinity, in practical salinity units
    :param pressure: Measured pressure, in decibars
    :param lon: Longitude
    :param lat: Latitude

    :return: Potential density in kg/m^3
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
    """Derive potential density from measured temperature, salinity, and
    pressure.

    :param temperature: Measure temperature, in degrees C
    :param conductivity: Measured conductivity, in mSiemens/cm
    :param pressure: Measured pressure, in decibars
    :param lon: Longitude
    :param lat: Latitude

    :return: Potential density in kg/m^3
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

    :param pressure: Numpy array of floats representing pressure, in
        dbar or psi
    :param latitude: Latitude (-90.0 to 90.0)
    :param depth_units: 'm' for meters, 'ft' for feet. Defaults to 'm'.
    :param pressure_units: 'dbar' for decibars, 'psi' for PSI. Defaults
        to 'dbar'.

    :return: A numpy array representing depth in meters or feet
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

    raw_oxygen_phase is expected to be in raw phase, raw_thermistor_temp
    in counts, pressure in dbar, and salinity in practical salinity (PSU)

    :param raw_oxygen_phase: SBE63 phase value, in microseconds
    :param thermistor_temp: SBE63 thermistor data to use are reference,
        in counts
    :param pressure: Converted pressure value from the attached CTD, in
        dbar
    :param salinity: Converted salinity value from the attached CTD, in
        practical salinity PSU
    :param coefs (Oxygen63Coefficients): calibration coefficients for
        the SBE63 sensor

    :return: converted Oxygen value, in ml/l
    """
    if thermistor_units == "volts":
        temperature = convert_sbe63_thermistor(thermistor, thermistor_coefs)
    elif thermistor_units == "C":
        temperature = thermistor
    else:
        raise ValueError

    oxygen_volts = raw_oxygen_phase / OXYGEN_PHASE_TO_VOLTS  # from the manual

    ksv = coefs.c0 + coefs.c1 * temperature + coefs.c2 * temperature**2

    # The following correction coefficients are all constants
    sol_b0 = -6.24523e-3
    sol_b1 = -7.37614e-3
    sol_b2 = -1.0341e-2
    sol_b3 = -8.17083e-3
    sol_c0 = -4.88682e-7

    ts = np.log((KELVIN_OFFSET_25C - temperature) / (KELVIN_OFFSET_0C + temperature))
    s_corr_exp = (
        salinity * (sol_b0 + sol_b1 * ts + sol_b2 * ts**2 + sol_b3 * ts**3) + sol_c0 * salinity**2
    )
    s_corr = e**s_corr_exp

    # temperature in Kelvin
    temperature_k = temperature + KELVIN_OFFSET_0C
    p_corr_exp = (coefs.e * pressure) / temperature_k
    p_corr = e**p_corr_exp

    # fmt: off
    ox_val = (
        (((coefs.a0 + coefs.a1 * temperature + coefs.a2 * oxygen_volts**2)
        / (coefs.b0 + coefs.b1 * oxygen_volts) - 1.0) / ksv) * s_corr * p_corr
    )
    # fmt: on

    return ox_val


def convert_sbe63_thermistor(
    instrument_output: np.ndarray,
    coefs: Thermistor63Coefficients,
):
    """Converts a SBE63 thermistor raw output array to temperature in
    ITS-90 deg C.

    :param instrument_output: raw values from the thermistor
    :param coefs: calibration coefficients for the thermistor in the
        SBE63 sensor

    :return: converted thermistor temperature values in ITS-90 deg C
    """
    log_raw = np.log((100000 * instrument_output) / (3.3 - instrument_output))
    temperature = (
        1 / (coefs.ta0 + coefs.ta1 * log_raw + coefs.ta2 * log_raw**2 + coefs.ta3 * log_raw**3)
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

    voltage is expected to be in volts, temperature in deg c, pressure
    in dbar, and salinity in practical salinity (PSU). All equation
    information comes from the June 2013 revision of the SBE43 manual

    :param voltage: SBE43 voltage
    :param temperature: temperature value converted to deg C
    :param pressure: Converted pressure value from the attached CTD, in
        dbar
    :param salinity: Converted salinity value from the attached CTD, in
        practical salinity PSU
    :param coefs: calibration coefficients for the SBE43 sensor
    :param apply_tau_correction: whether or not to run tau correction
    :param apply_hysteresis_correction: whether or not to run hysteresis
        correction
    :param window_size: size of the window to use for tau correction, if
        applicable, in seconds
    :param sample_interval: sample rate of the data to be used for tau
        correction, if applicable. In seconds.

    :return: converted Oxygen values, in ml/l
    """
    # start with all 0 for the dvdt
    dvdt_values = np.zeros(len(voltage))
    if apply_tau_correction:
        # Calculates how many scans to have on either side of our median
        # point, accounting for going out of index bounds
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

    voltage is expected to be in volts, temperature in deg c, pressure
    in dbar, and salinity in practical salinity (PSU). All equation
    information comes from the June 2013 revision of the SBE43 manual.
    Expects that hysteresis correction is already performed on the
    incoming voltage, if desired.

    :param voltage: SBE43 voltage
    :param temperature: temperature value converted to deg C
    :param pressure: Converted pressure value from the attached CTD, in
        dbar
    :param salinity: Converted salinity value from the attached CTD, in
        practical salinity PSU
    :param coefs: calibration coefficients for the SBE43 sensor
    :param dvdt_value: derivative value of voltage with respect to time
        at this point. Expected to be 0 if not using Tau correction

    :return: converted Oxygen value, in ml/l
    """

    # Oxygen Solubility equation constants, From SBE43 Manual Appendix A
    a0 = 2.00907
    a1 = 3.22014
    a2 = 4.0501
    a3 = 4.94457
    a4 = -0.256847
    a5 = 3.88767
    b0 = -0.00624523
    b1 = -0.00737614
    b2 = -0.010341
    b3 = -0.00817083
    c0 = -0.000000488682

    ts = np.log((KELVIN_OFFSET_25C - temperature) / (KELVIN_OFFSET_0C + temperature))
    a_term = a0 + a1 * ts + a2 * ts**2 + a3 * ts**3 + a4 * ts**4 + a5 * ts**5
    b_term = salinity * (b0 + b1 * ts + b2 * ts**2 + b3 * ts**3)
    c_term = c0 * salinity**2
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

    From Application Note 64.

    :param ox_values: oxygen values, already converted to ml/L

    :return: oxygen values converted to milligrams/Liter
    """

    return ox_values * OXYGEN_MLPERL_TO_MGPERL


def convert_oxygen_to_umol_per_kg(ox_values: np.ndarray, potential_density: np.ndarray):
    """Converts given oxygen values to milligrams/kg.

    Note: Sigma-Theta is expected to be calculated via gsw_sigma0,
    meaning is it technically potential density anomaly. Calculating
    using gsw_rho(SA, CT, p_ref = 0) results in actual potential
    density, but this function already does the converison, so values
    will need to have 1000 subtracted from them before being passed into
    this function. The function is done this way to stay matching to
    Application Note 64, but the results of either method are identical.

    :param ox_values: oxygen values, already converted to ml/L
    :param potential_density: potential density (sigma-theta) values.
        Expected to be the same length as ox_values

    :return: oxygen values converted to milligrams/Liter
    """

    oxygen_umolkg = (ox_values * OXYGEN_MLPERL_TO_UMOLPERKG) / (potential_density + 1000)
    return oxygen_umolkg


def convert_eco(
    raw: np.ndarray,
    coefs: ECOCoefficients,
):
    """Converts a raw value for any ECO measurand.

    :param raw: raw counts for digital, raw volts for analog
    :param coefs (ChlorophyllACoefficients): calibration coefficients

    :return: converted ECO measurement in calibration units
    """
    converted = coefs.slope * (raw - coefs.offset)

    return converted


def convert_sbe18_ph(
    raw_ph: np.ndarray,
    temperature: np.ndarray,
    coefs: PH18Coefficients,
):
    """Converts a raw voltage value for pH.

    All equation information comes from application note 18-1

    :param raw_ph: raw output voltage from pH sensor (0-5V)
    :param temperature: temperature value to use for temperature
        compensation in degrees C
    :param coefs: slope and offset for the pH sensor

    :return: converted pH
    """
    ph = 7 + (raw_ph - coefs.offset) / (
        1.98416e-4 * (temperature + KELVIN_OFFSET_0C) * coefs.slope
    )
    return ph


def convert_par_logarithmic(
    raw_par: np.ndarray,
    coefs: PARCoefficients,
):
    """Converts a raw voltage value for PAR to µmol photons/m2*s.

    All equation information comes from application note 96

    :param raw_par: raw output voltage from PAR sensor
    :param coefs: calibration coefficients for the PAR sensor

    :return: converted PAR in µmol photons/m2*s
    """
    par = coefs.multiplier * coefs.im * 10 ** ((raw_par - coefs.a0) / coefs.a1)

    return par


def convert_nitrate(
    volts: np.ndarray, dac_min: float, dac_max: float, units: Literal["uMNO3", "mgNL"] = "uMNO3"
):
    """Convert SUNA raw voltages to uMNO3 or mgNL

    :param volts: raw output voltage from a SUNA
    :param dac_min: NO3 value that corresponds to v_min
    :param dac_max: NO3 value that corresponds to v_max
    :param units: conversion output units, defaults to 'uMNO3'
    :return: converted nitrate
    """
    v_min = 0.095
    v_max = 4.095
    a1 = (dac_min - dac_max) / (v_min - v_max)
    a0 = dac_max - v_max * a1

    nitrate = a1 * volts + a0

    if units == "mgNL":
        nitrate *= UMNO3_TO_MGNL

    return nitrate


def convert_ph_voltage_counts(ph_counts: np.ndarray):
    """Convert pH voltage counts to a floating point value

    :param ph_counts: pH voltage counts
    :return: pH voltage
    """
    adc_vref = 2.5
    gain = 1
    adc_23bit = 8388608
    ph_volts = adc_vref / gain * (ph_counts / adc_23bit - 1)
    return ph_volts


def _calculate_nernst(temperature: np.ndarray):
    """Calculate the nernst term using natual log

    :param temperature: temperature in kelvin
    :return: the nernst term (J/Coulomb; electrical potential; volts)
    """
    nernst_term = R * temperature / F * np.log(10)
    return nernst_term


def convert_internal_seafet_ph(
    ph_counts: np.ndarray,
    temperature: np.ndarray,
    coefs: PHSeaFETInternalCoefficients,
):
    """Calculates the internal pH on the total scale given the
    temperature and internal FET voltage

    :param ph_counts: pH voltage counts
    :param temperature: sample temperature
    :param coefs: SeaFET calibration coefficients
    :return: calculated pH on the total scale for the SeaFET internal
        reference
    """
    ph_volts = convert_ph_voltage_counts(ph_counts)

    # Eo(T) or temperature offset
    temperature_offset = coefs.k2 * temperature

    # Eo the cell reference voltage at in-situ conditions
    cell_ref_volts = coefs.k0 + temperature_offset
    nernst_term = _calculate_nernst(temperature + KELVIN_OFFSET_0C)
    ph = (ph_volts - cell_ref_volts) / nernst_term
    return ph


# pylint: disable=too-many-statements #TODO: break this function up
def convert_external_seafet_ph(
    ph_counts: np.ndarray,
    temperature: np.ndarray,
    salinity: np.ndarray,
    pressure: np.ndarray,
    coefs: PHSeaFETExternalCoefficients,
):
    """Calculates the external pH on the total scale given temperature,
    salinity, pressure and FET voltage counts

    :param ph_counts: ISFET external voltage counts
    :param temperature: sample temperature in Celsius
    :param salinity: sample salinity in psu
    :param pressure: sample pressure in dbar
    """

    ts_cr = 0.14  # relative concentration of sulfate in SW
    ts_mm = 96.062  # [g/mol] molar mass of sulfate
    cl_cr = 0.99889  # relative concentration of chloride in SW
    cl_mm = 35.453  # [g/mol] molar mass of chloride
    cl_to_s = 1.80655  # [ppt, 10^{-3}] Chlorinity to Salinity

    def _molar_concentration(concentration: float, molar_mass: float, salinity: np.ndarray):
        """
        An estimate of constituent concentration in seawater is made
        from salinity on the basis of contancy of composition.

        :param concentration: Relative concentration of constituent
        :param molar_mass: Molar mass of constituent
        :param salinity: Salinity in PSU
        :return: molar concentration
        """
        return (concentration / molar_mass) * (salinity / cl_to_s)

    def _molar_conc_chloride(salinity: np.ndarray):
        """Calculates the molar concentration of chloride

        :param salinity: Salinity in PSU
        :return: molar concentration of chloride
        """
        return _molar_concentration(cl_cr, cl_mm, salinity) * 1000 / (1000 - 1.005 * salinity)

    def _molar_conc_sulfate(salinity: np.ndarray):
        """Calculates the molar concentration of total sulfate

        :param salinity: Salinity in PSU
        :return: molar concentration of sulfate
        """
        return _molar_concentration(ts_cr, ts_mm, salinity)

    def _calculate_ionic_strength(salinity: np.ndarray):
        """Compute Salinity Ionic strength
        Ionic Strength (mol/kg H2O) from Dickson "Guide to Best
        Practices for Ocean CO2 Measurements", 2007, Chapter 5, page 11

        :param salinity: Salinity in PSU
        :return: Salinity ionic strength
        """
        c0 = 19.924
        c1 = 1000
        c2 = 1.005

        # 19.924*(S) / (1000 - 1.005*(S))
        ionic_strength = c0 * salinity / (c1 - c2 * salinity)
        return ionic_strength

    def _calculate_ks(
        temperature: np.ndarray,
        ionic_strength: np.ndarray,
        salinity: np.ndarray,
        pressure: np.ndarray,
    ):
        """Dissociation constant of sulfuric acid in seawater
        Dickson, A. G., J. Chemical Thermodynamics, 22:113-127, 1990
        The goodness of fit is .021.
        It was given in mol/kg-H2O. I convert it to mol/kg-SW.
        TYPO on p. 121: the constant e9 should be e8.
        This is from eqs 22 and 23 on p. 123, and Table 4 on p 121:

        :param temperature: Temperature in kelvin
        :param ionic_strength: Ionic strength
        :param salinity: Salinity in PSU
        :param pressure: Pressure in dbar
        :return: Dissociation constant of sulfuric acid in seawater
        """

        # *********** this should be re-tested
        c0 = -4276.1
        c1 = 141.328
        c2 = -23.093
        c3 = -13856
        c4 = 324.57
        c5 = -47.986
        c6 = 35474
        c7 = -771.54
        c8 = 114.723
        c9 = -2698
        c10 = 1776

        ln_ks = (
            c0 / temperature
            + c1
            + c2 * np.log(temperature)
            + (c3 / temperature + c4 + c5 * np.log(temperature)) * np.sqrt(ionic_strength)
            + (c6 / temperature + c7 + c8 * np.log(temperature)) * ionic_strength
            + (c9 / temperature) * np.sqrt(ionic_strength) * ionic_strength
            + (c10 / temperature) * np.pow(ionic_strength, 2)
        )

        # this is on the free pH scale in mol/kg-H2O
        khso4 = np.exp(ln_ks) * (1 - 0.001005 * salinity)  # convert to mol/kg-SW

        # UCI has two calculateKS functions. The first returns khso4 here,
        # but in practice the second is always used which adds the following

        temperature_c = temperature - KELVIN_OFFSET_0C
        # Partial molal volume and compressibility change for HSO4
        delta_vhso4 = -18.03 + 0.0466 * temperature_c + 0.000316 * temperature_c**2
        kappa_hso4 = (-4.53 + 0.09 * temperature_c) / 1000

        #  per Yui Press changed from dbar to bar here by / 10
        ln_khso4_fac = (
            (-delta_vhso4 + 0.5 * kappa_hso4 * (pressure / 10))
            * (pressure / 10)
            / (R * 10 * temperature)
        )

        #  bisulfate association constant at T, S, P
        khso4_tps = khso4 * np.exp(ln_khso4_fac)

        return khso4_tps

    def _calculate_adh(temperature):
        """Calculates the Debeye-Huckel constant (temperature [Celcius]
        dependence only)

        :param temperature: temperature in C
        :return: Debeye-Huckel constant
        """
        # This fit was made by Ken Johnson from data presented by:
        # Khoo et al. (Anal. Chem., 49, 29-34, 1977).
        c0 = 0.49172143
        c1 = 0.00067524
        # Modified from 0.00000343 to 0.0000034286,
        # email from Ken with newest version of MBARI code by Charles Branham 3/9/16
        c2 = 0.0000034286

        # 0.00000343*(t)*(t) + 0.00067524*(t) + 0.49172143
        adh = c0 + c1 * temperature + c2 * temperature**2

        return adh

    def _calculate_log_gamma_hcl(
        adh: np.ndarray, ionic_strength: np.ndarray, temperature: np.ndarray
    ):
        # pylint: disable=anomalous-backslash-in-string
        """
        Khoo et al. (Anal. Chem., 49, 29-34, 1977).
        log γ±(HCl) = -A·√I / (1 + ρ·√I) + (B₀ + B₁·T)·I
        As implemented in the calibration .xls

        :param adh: Debeye-Huckel constant
        :param ionic_strength: ionic strength
        :param temperature: Temperatur in C
        :return: Log Gamma HCL
        """
        # NOTE: There is a disagreement between Ken Johnson and Yui Takeshita
        # as to whether in-situ pressure correction is needed or already taken care of
        # by the ThermPress term being added to E0. In email correspondence
        # with Dave Murphy, Ken Johnson finally stated that no pressure correction
        # should be applied to Gamma_HCl

        rho = 1.394
        b0 = 0.08885
        b1 = 0.000111

        # As per instructions of KJ, should consider the data in
        # Dickson (J. Chem. Thermodynamics, 22, 113-127, 1990)
        log_gamma_hcl = (
            -1 * adh * np.sqrt(ionic_strength) / (1 + rho * np.sqrt(ionic_strength))
            + (b0 - b1 * temperature) * ionic_strength
        )

        return log_gamma_hcl

    def _calculate_thermal_pressure(temperature: np.ndarray, pressure: np.ndarray):
        """ThermPress at in-situ
        Thermpress = (-V_Cl x P + 0.5 K_Cl x P^2)/10F
        Where,
        P: Pressure in 'bar'
        V_Cl: Chloride partial molal volume (in cm^3 mol^-1)
        K_Cl: Chloride partial molal compressibility
        (in cm^3 mol^-1 bar^-1) (can be neglected [4, pg. 876])

        :param temperature: temperature in C
        :param pressure: pressure in dbar
        :return: thermal pressure
        """
        # // Partial molal volume and compressibility change for HCl from Millero
        delta_vhcl = 17.85 + 0.1044 * temperature - 0.001316 * temperature**2

        # // Pressure is in 'dbar' and need to divide by 10 to convert to 'bar'
        pressure_bar = pressure / 10

        # // Thermpress term
        thermal_pressure = -delta_vhcl * 0.0242 / (23061 * 1.01) * pressure_bar

        return thermal_pressure

    # Intermediate terms
    # tempK;   # Sample Temperature (K)
    temperature_k = temperature + KELVIN_OFFSET_0C
    # ST;        # Nernst Term
    st = _calculate_nernst(temperature=temperature_k)
    # mCl;       # Molar concentration of chloride
    mcl = _molar_conc_chloride(salinity=salinity)
    # TSO4;      # Molar concentration of total sulfate
    tso4 = _molar_conc_sulfate(salinity=salinity)
    # IonS;      # Salinity Ionic strength
    ionic_strength = _calculate_ionic_strength(salinity=salinity)

    adh = _calculate_adh(temperature=temperature)

    log_gamma_hcl = _calculate_log_gamma_hcl(
        adh=adh, ionic_strength=ionic_strength, temperature=temperature
    )

    thermal_pressure = _calculate_thermal_pressure(temperature=temperature, pressure=pressure)

    # // Dissociation constant of sulfuric acid in seawater
    # KSO4;      # Dissociation constant of sulfuric acid in seawater
    # this.KSO4 = calculateKS(this.tempK, this.IonS, salinity, tempC, pressure);
    kso4 = _calculate_ks(
        temperature=temperature_k,
        ionic_strength=ionic_strength,
        salinity=salinity,
        pressure=pressure,
    )

    ph_volts = convert_ph_voltage_counts(ph_counts)

    # Eo(T) or temperature offset
    eot = coefs.k2 * temperature

    # Eo(P) or pressure offset
    eop = (
        coefs.f1 * pressure
        + coefs.f2 * pressure**2
        + coefs.f3 * pressure**3
        + coefs.f4 * pressure**4
        + coefs.f5 * pressure**5
        + coefs.f6 * pressure**6
    )

    # Calculate the External pH in Free Scale, firmware applies 2 factor here.
    # Convert pHfree(mol/kg H2O) to pHfree(mol/kg sln) per KSJ by using the
    # ratio (1000 - 1.005 * Salinity)/1000
    free_scale_ph = (
        (ph_volts - eot - eop - thermal_pressure - coefs.k0) / st
        + np.log10(mcl)
        + (2 * log_gamma_hcl)
        - np.log10((1000 - 1.005 * salinity) / 1000)
    )

    total_ph = free_scale_ph - np.log10(1 + tso4 / kso4)

    return total_ph


def convert_internal_seafet_temperature(temperature_counts: np.ndarray):
    """Converts the raw internal temperature counts to degrees Celcius

    :param temperature_counts: raw internal temperature counts
    :return: internal temperature in Celcius
    """
    slope = 175.72
    offset = -46.85
    int_16bit = 2**16
    temperature = temperature_counts / int_16bit * slope + offset

    return temperature


def convert_seafet_relative_humidity(humidity_counts: np.ndarray, temperature: np.ndarray):
    """Convert relative humidity counts to percent

    :param humidity_counts: raw relative humidity counts
    :param temperature: converted internal temperature in Celcius
    :return: temperature compensated relative humidity in percent
    """
    slope = 125
    offset = -6
    int_16bit = 2**16
    max_humidity = 119
    temperature_coefficient = -0.15
    temperature_25c = 25

    # Uncompensated relative humidity
    relative_humidity = slope * humidity_counts / int_16bit + offset

    for n, humidity in enumerate(relative_humidity):
        # Theoretically, uncompensated relative humidity can be up to 119%
        if 0 <= humidity < max_humidity:
            relative_humidity[n] = humidity + temperature_coefficient * (
                temperature_25c - temperature[n]
            )

    np.clip(relative_humidity, a_min=0, a_max=100)

    return relative_humidity
