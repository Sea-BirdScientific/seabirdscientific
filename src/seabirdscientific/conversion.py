#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""A collection of raw data conversion functions."""



# Native imports
from math import e, floor
from typing import Literal
import warnings

# Third-party imports
import gsw
import numpy as np
from numpy.polynomial import Polynomial
from scipy import stats

# Sea-Bird imports

# Internal imports
import seabirdscientific.cal_coefficients as cc


DBAR_TO_PSI = 1.450377
PSI_TO_DBAR = 0.6894759
OXYGEN_PHASE_TO_VOLTS = 39.457071
KELVIN_OFFSET_0C = 273.15
KELVIN_OFFSET_25C = 298.15
OXYGEN_MLPERL_TO_MGPERL = 1.42903
OXYGEN_MLPERL_TO_UMOLPERKG = 44660
# taken from https://blog.seabird.com/ufaqs/what-is-the-difference-in-temperature-expressions-between-ipts-68-and-its-90/
ITS90_TO_IPTS68 = 1.00024
# micro moles of nitrate to milligrams of nitrogen per liter
UMNO3_TO_MGNL = 0.014007
# [J K^{-1} mol^{-1}] Gas constant from SBS application note 99
R = 8.3144621
# [Coulombs mol^{-1}] Faraday constant from SBS application note 99
F = 96485.365


def convert_temperature(
    temperature_counts_in: np.ndarray,
    coefs: cc.TemperatureCoefficients,
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

    :return: temperature val converted to ITS-90 or IPTS68 in degrees C or F
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


def convert_temperature_frequency(
    frequency: np.ndarray,
    coefs: cc.TemperatureFrequencyCoefficients,
    standard: Literal["ITS90", "IPTS68"] = "ITS90",
    units: Literal["C", "F"] = "C",
):
    """Convert raw frequency to temperature in degrees Celsius or degrees Fahrenheit

    :param frequency: raw frequency from the temperature sensor
    :param coefs: calibration coefficients for the temperature sensor
    :return: temperature in Celsius or Fahrenheit
    """
    fLog = np.log(coefs.f0 / frequency)
    temperature = (
        1 / (coefs.g + coefs.h * fLog + coefs.i * fLog**2 + coefs.j * fLog**3) - KELVIN_OFFSET_0C
    )

    if standard == "IPTS68":
        temperature *= ITS90_TO_IPTS68
    if units == "F":
        temperature = temperature * 9 / 5 + 32  # Convert C to F

    return temperature


def convert_pressure(
    pressure_count: np.ndarray,
    compensation_voltage: np.ndarray,
    coefs: cc.PressureCoefficients,
    units: Literal["dbar", "psia", "psig"] = "psig",
):
    """Converts pressure counts to sea pressure (psig and dbar) and absolute pressure (psia)

    pressure_count and compensation_voltage are expected to be raw data
    from an instrument in A/D counts

    :param pressure_count: pressure value to convert, in A/D counts
    :param compensation_voltage: pressure temperature compensation
        voltage, in counts or volts depending on the instrument
    :param coefs: calibration coefficients for the pressure sensor
    :param units: whether or not to use psig or dbar as the returned
        unit type

    :return: sea pressure val in dbar or PSIG
    """
    sea_level_pressure = 14.7

    t = (
        coefs.ptempa0
        + coefs.ptempa1 * compensation_voltage
        + coefs.ptempa2 * compensation_voltage**2
    )
    x = pressure_count - coefs.ptca0 - coefs.ptca1 * t - coefs.ptca2 * t**2
    n = x * coefs.ptcb0 / (coefs.ptcb0 + coefs.ptcb1 * t + coefs.ptcb2 * t**2)
    pressure = coefs.pa0 + coefs.pa1 * n + coefs.pa2 * n**2

    if units == "dbar" or units == "psig":
        pressure -= sea_level_pressure

    if units == "dbar":
        pressure *= PSI_TO_DBAR

    return pressure


def convert_pressure_digiquartz(
    pressure_count: np.ndarray,
    compensation_voltage: np.ndarray,
    coefs: cc.PressureDigiquartzCoefficients,
    units: Literal["dbar", "psia"],
    sample_interval: float,
):
    """Converts pressure counts to PSIA (pounds per square inch, abolute) or dbar for a digiquartz pressure sensor.

    pressure_count and compensation_voltage are expected to be raw data
    from an instrument in A/D counts

    :param pressure_count: pressure value to convert, in A/D counts
    :param compensation_voltage: pressure temperature compensation
        voltage, in counts or volts depending on the instrument
    :param coefs: calibration coefficients for the digiquartz pressure sensor
    :param units: whether or not to use psia or dbar as the returned
        unit type
    :param sample_interval: sample rate of the data to be used for temperature compensation correction, in seconds
    :return: pressure val in PSIA or dbar
    """
    sea_level_pressure = 14.7
    # First, average temperature compensation over 30 seconds
    max_scans_in_30_seconds = 720
    scans_in_window = floor(30 / sample_interval)
    scans_in_window = max(scans_in_window, 1)
    scans_in_window = min(scans_in_window, max_scans_in_30_seconds)

    rolling_sum = compensation_voltage[0] * scans_in_window
    modified_compensation_voltage = compensation_voltage.copy()

    for i in range(0, len(compensation_voltage)):
        if i < scans_in_window:
            # remove a copy of 0-index value from rolling sum
            rolling_sum -= compensation_voltage[0]
        else:
            # remove oldest value from rolling sum
            rolling_sum -= compensation_voltage[i - scans_in_window]

        rolling_sum += compensation_voltage[i]
        modified_compensation_voltage[i] = (
            rolling_sum / scans_in_window * coefs.AD590M + coefs.AD590B
        )

    # Now, calculate pressure

    t = 1 / pressure_count * 1000000  # convert to period in usec
    c = (
        coefs.c1
        + coefs.c2 * modified_compensation_voltage
        + coefs.c3 * modified_compensation_voltage**2
    )
    d = coefs.d1 + coefs.d2 * modified_compensation_voltage
    t0 = (
        coefs.t1
        + coefs.t2 * modified_compensation_voltage
        + coefs.t3 * modified_compensation_voltage**2
        + coefs.t4 * modified_compensation_voltage**3
        + coefs.t5 * modified_compensation_voltage**4
    )

    t0_squared_over_t_squared = (t0**2) / (t**2)
    one_minus_ratio = 1 - t0_squared_over_t_squared
    p = c * one_minus_ratio * (1 - d * one_minus_ratio)
    abs_pressure = p - sea_level_pressure
    if units == "dbar":
        abs_pressure *= PSI_TO_DBAR
    return abs_pressure


def convert_conductivity(
    conductivity_count: np.ndarray,
    temperature: np.ndarray,
    pressure: np.ndarray,
    coefs: cc.ConductivityCoefficients,
    scalar: float = 1.0,
):
    """Converts raw conductivity counts to S/m.

    Data is expected to be raw data from instrument in A/D counts

    :param conductivity_count: conductivity value to convert, in A/D
        counts
    :param temperature: reference temperature, in degrees C
    :param pressure: reference pressure, in dbar
    :param coefs: calibration coefficient for the conductivity sensor
    :param scalar: value to multiply by at the end. For most instruments, this is 1. For SBE911, it is 1/10

    :return: conductivity val converted to S/m
    """
    f = conductivity_count * np.sqrt(1 + coefs.wbotc * temperature) / 1000
    numerator = coefs.g + coefs.h * f**2 + coefs.i * f**3 + coefs.j * f**4
    denominator = 1 + coefs.ctcor * temperature + coefs.cpcor * pressure
    return numerator / denominator * scalar


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
    coefs: cc.Oxygen63Coefficients,
    thermistor_coefs: cc.Thermistor63Coefficients,
    thermistor_units: Literal["volts", "C"] = "volts",  # Is this volts or frequency?
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
    :param coefs (cc.Oxygen63Coefficients): calibration coefficients for
        the SBE63 sensor
    :param thermistor_coefs (cc.Thermistor63Coefficients): calibration coefficients for
        the SBE63 thermistor sensor
    :param thermistor_units: units of thermistor_temp input

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
    coefs: cc.Thermistor63Coefficients,
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
    coefs: cc.Oxygen43Coefficients,
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
    coefs: cc.Oxygen43Coefficients,
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
    coefs: cc.ECOCoefficients,
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
    coefs: cc.PH18Coefficients,
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
    volts: np.ndarray,
    coefs: cc.PARCoefficients,
):
    """Converts a raw voltage value for underwater PAR.

    All equation information comes from application note 96

    conversion_factor = 1.0 for units of μmol photons/m2*s

    :param raw_par: raw output voltage from PAR sensor
    :param coefs: calibration coefficients for the PAR sensor

    :return: converted PAR in µmol photons/m2*s
    """
    exponent = (volts - coefs.a0) / coefs.a1
    par = coefs.multiplier * coefs.im * 10**exponent

    return par


def convert_spar_logarithmic(
    volts: np.ndarray,
    coefs: cc.SPARCoefficients,
):
    """Converts a raw voltage value for logarithmic surface PAR.

    All equation information comes from application note 96

    conversion_factor = 1.0 for units of μmol photons/m2*s

    :param volts: raw output voltage from SPAR sensor
    :param coefs: coefficients for the SPAR sensors

    :return: converted surface PAR in µmol photons/m2*s
    """
    exponent = (volts - coefs.a0) / coefs.a1
    spar = coefs.conversion_factor * coefs.im * 10**exponent

    return spar


def convert_spar_linear(
    volts: np.ndarray,
    coefs: cc.SPARCoefficients,
):
    """Converts a raw voltage value for linear surface PAR.

    All equation information comes from application note 96

    conversion_factor = 1.0 for units of μmol photons/m2*s

    :param volts: raw output voltage from SPAR sensor
    :param coefs: coefficients for the SPAR sensors

    :return: converted surface PAR in µmol photons/m2*s
    """
    spar = coefs.im * coefs.a1 * (volts - coefs.a0) * coefs.conversion_factor

    return spar


def convert_spar_biospherical(
    volts: np.ndarray,
    coefs: cc.SPARCoefficients,
):
    """Converts a raw voltage value for biospherical surface PAR.

    All equation information comes from application note 11S

    :param volts: raw output voltage from SPAR sensor
    :param coefs: coefficients for the SPAR sensors

    :return: converted surface PAR in µmol photons/m2*s
    """
    spar = volts * coefs.conversion_factor

    return spar


def convert_nitrate(
    volts: np.ndarray,
    dac_min: float,
    dac_max: float,
    units: Literal["uMNO3", "mgNL"] = "uMNO3",
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


def _calculate_nernst(temperature: np.ndarray) -> np.ndarray:
    """Calculate the nernst term using natual log

    :param temperature: temperature in kelvin
    :return: the nernst term (J/Coulomb; electrical potential; volts)
    """
    nernst_term = R * temperature * np.log(10) / F
    return nernst_term


def convert_internal_seafet_ph(
    raw_ph: np.ndarray = 0,
    temperature: np.ndarray = 0,
    coefs: cc.PHSeaFETInternalCoefficients = cc.PHSeaFETInternalCoefficients(),
    ph_units: Literal["counts", "volts"] = "counts",
    ph_counts: np.ndarray = None,
):
    """Calculates the internal pH on the total scale given the
    temperature and internal FET voltage

    :param raw_ph: Raw voltage or voltage counts
    :param temperature: Sample temperature
    :param coefs: SeaFET calibration coefficients
    :param ph_units: The units of raw_ph, defaults to 'counts'
    :param ph_counts: Deprecated. pH voltage counts
    :return: calculated pH on the total scale for the SeaFET internal
        reference
    """
    if ph_counts is not None:
        warnings.warn("Deprecated, use raw_ph", DeprecationWarning)
        ph_volts = convert_ph_voltage_counts(ph_counts)
    elif ph_units == 'counts':
        ph_volts = convert_ph_voltage_counts(raw_ph)
    else:  # ph_counts == 'volts'
        ph_volts = raw_ph

    nernst = _calculate_nernst(temperature + KELVIN_OFFSET_0C)
    ph = (ph_volts - coefs.kdf0 - coefs.kdf2 * temperature) / nernst
    return ph


def _calculate_thermal_pressure(temperature: np.ndarray, pressure: np.ndarray):
    """ThermPress at in-situ
    Thermpress = (-V_Cl x P + 0.5 K_Cl x P^2)/10F
    Where,
    P: Pressure in 'bar'
    V_Cl: Chloride partial molal volume (in cm^3 mol^-1)
    K_Cl: Chloride partial molal compressibility
    (in cm^3 mol^-1 bar^-1) (can be neglected [4, pg. 876])

    :param temperature: temperature in C
    :param pressure: pressure in bar
    :return: thermal pressure
    """
    # // Partial molal volume and compressibility change for HCl from Millero
    delta_vhcl = _partial_molal_hcl_volume(temperature)

    # // Thermpress term
    thermal_pressure = -delta_vhcl * 0.0242 / (23061 * 1.01) * pressure

    return thermal_pressure


def _total_chloride_in_seawater(salinity: np.ndarray) -> np.ndarray:
    """From SBS application note 99. Calculated as (Dickson et al. 2007)

    :param salinity: Salinity in PSU
    :return: Total chloride in seawater
    """
    # 0.99889  relative concentration of chloride in SW
    # 35.453  [g/mol] molar mass of chloride
    # 1.80655  [ppt, 10^{-3}] Chlorinity to Salinity
    factor_1 = 0.99889 / 35.453
    factor_2 = salinity / 1.80655
    factor_3 = 1 / (1 - 1.005e-3 * salinity)
    total_chloride = factor_1 * factor_2 * factor_3
    return total_chloride


def _total_sulfate_in_seawater(salinity: np.ndarray):
    """From SBS application note 99. Calculated as (Dickson et al. 2007)

    :param salinity: Salinity in PSU
    :return: Total sulfate in seawater
    """
    # 0.14  relative concentration of sulfate in SW
    # 96.062  [g/mol] molar mass of sulfate
    # 1.80655  [ppt, 10^{-3}] Chlorinity to Salinity
    total_sulfate = (0.14 / 96.062) * (salinity / 1.80655)
    return total_sulfate


def _sample_ionic_strength(salinity: np.ndarray) -> np.ndarray:
    """From SBS application note 99. The sample ionic strength is
    calculated as (Dickson et al. 2007)

    :param salinity: Salinity in PSU
    :return: Sample ionic strength
    """
    ionic_strength = (19.924 * salinity) / (1000 - 1.005 * salinity)
    return ionic_strength


def _debye_huckel_constant_for_hcl_activity(temperature: np.ndarray):
    """From SBS application note 99. This constant is calculated as
    (Khoo et al. 1977)

    :param temperature: Temperature in degrees C
    :return: Debye-Huckel constant for activity of HCl
    """
    activity = 3.4286e-6 * temperature**2 + 6.7524e-4 * temperature + 0.49172143
    return activity


def _log_of_hcl_activity_coefficient_of_t(salinity: np.ndarray, temperature: np.ndarray):
    """From SBS application note 99. Calculated as (Khoo et al. 1977)

    :param salinity: Salinity in PSU
    :param temperature: Temperature in degrees C
    :return: Logarithm of HCl activity coefficient as a function of temperature
    """
    i = _sample_ionic_strength(salinity)
    a_dh = _debye_huckel_constant_for_hcl_activity(temperature)
    term_1 = (-a_dh * np.sqrt(i)) / (1 + 1.394 * np.sqrt(i))
    term_2 = (0.08885 - 1.11e-4 * temperature) * i
    log_y_hcl = term_1 + term_2
    return log_y_hcl


def _log_of_hcl_activity_coefficient_of_tp(
    salinity: np.ndarray,
    temperature: np.ndarray,
    pressure: np.ndarray,
):
    """From SBS application note 99. Calculated as (Johnson et al. 2017)

    :param salinity: Salinity in PSU
    :param temperature: Temperature in degrees C
    :param pressure: Pressure in bar
    :return: Logarithm of HCl activity coefficient as a function of
        temperature and pressure
    """
    log_y_hcl = _log_of_hcl_activity_coefficient_of_t(salinity, temperature)
    v_hcl = _partial_molal_hcl_volume(temperature)
    t_kelvin = temperature + KELVIN_OFFSET_0C
    term_2 = (v_hcl * pressure) / (np.log(10) * R * t_kelvin * 10) / 2
    log_y_hcl_tp = log_y_hcl + term_2
    return log_y_hcl_tp


def _acid_dissociation_constant_of_hso4(salinity: np.ndarray, temperature: np.ndarray):
    """From SBS application note 99. Calculated as (Dickson et al. 2007)

    :param salinity: Salinty in PSU
    :param temperature: Temperature in Kelvin
    :return: _description_
    """
    i = _sample_ionic_strength(salinity)
    term_1 = -4276.1 / temperature + 141.328 - 23.093 * np.log(temperature)
    term_2 = (-13856 / temperature + 324.57 - 47.986 * np.log(temperature)) * np.sqrt(i)
    term_3 = (35474 / temperature - 771.54 + 114.723 * np.log(temperature)) * i
    term_4 = (2698 / temperature) * i**1.5
    term_5 = (1776 / temperature) * i**2
    k_s = (1 - 1.005e-3 * salinity) * np.exp(term_1 + term_2 + term_3 - term_4 + term_5)
    return k_s


def _acid_dissociation_constant_of_hso4_tp(
    salinity: np.ndarray,
    temperature: np.ndarray,
    pressure: np.ndarray,
):
    """From SBS application note 99. Calculated as (Millero 1982)

    :param salinity: Salinity in PSU
    :param temperature: Temperature in degrees C
    :param pressure: Pressure in bar
    :return: Acid dissociation constant of HSO4
    """
    t_kelvin = temperature + KELVIN_OFFSET_0C
    k_s = _acid_dissociation_constant_of_hso4(salinity, t_kelvin)
    v_bar_s = _partial_molal_hso4_volume(temperature)
    k_bar_s = _hso4_compressibility(temperature)
    exponent = (-v_bar_s * pressure + 0.5 * k_bar_s * pressure**2) / (R * t_kelvin * 10)
    k_stp = k_s * np.exp(exponent)
    return k_stp


def _pressure_response(pressure: np.ndarray, coefs: cc.PHSeaFETExternalCoefficients):
    """The sensor pressure response function from SBS application note
    99.

    :param pressure: Pressure in bar
    :param coefs: External pH coefficients
    :return: The pressure response
    """
    term_1 = coefs.f1 * pressure
    term_2 = coefs.f2 * pressure**2
    term_3 = coefs.f3 * pressure**3
    term_4 = coefs.f4 * pressure**4
    term_5 = coefs.f5 * pressure**5
    term_6 = coefs.f6 * pressure**6
    return term_1 + term_2 + term_3 + term_4 + term_5 + term_6


def _partial_molal_hcl_volume(temperature: np.ndarray):
    """From SBS application note 99. Calculated as (Millero 1983)
    Note: AN99 has a typo and is missing temperature in the second term

    :param temperature: Temperature in degrees C
    :return: Partial Molal Volume of HCl
    """
    volume = 17.85 + 0.1044 * temperature - 1.316e-3 * temperature**2
    return volume


def _partial_molal_hso4_volume(temperature: np.ndarray) -> np.ndarray:
    """From SBS application note 99. Calculated as (Millero 1983)

    :param temperature: Temperature in dgrees C
    :return: Partial Molal Volume of HSO4
    """
    volume = -18.03 + 0.0466 * temperature + 3.16e-4 * temperature**2
    return volume


def _hso4_compressibility(temperature: np.ndarray):
    """From SBS application note 99. Calculated as (Millero 1983)

    :param temperature: Temperature in degrees C
    :return: Compressibility of HSO4
    """
    compressibility = (-4.53 + 0.09 * temperature) / 1000
    return compressibility


def convert_external_seafet_ph(
    raw_ph: np.ndarray = 0,
    temperature: np.ndarray = 0,
    salinity: np.ndarray = 0,
    pressure: np.ndarray = 0,
    coefs: cc.PHSeaFETExternalCoefficients = cc.PHSeaFETExternalCoefficients(),
    ph_units: Literal["counts", "volts"] = "counts",
    formula_version: Literal['legacy', '1.3'] = '1.3',
    ph_counts: np.ndarray = None,
):
    """External pH for the SeaFET, SeapHOx, and Float. From SBS
    Application Note 99 and "Processing BGC-Argo pH data at the DAC
    level"
    
    https://www.seabird.com/asset-get.download.jsa?id=69833850609
    https://archimer.ifremer.fr/doc/00460/57195/

    :param raw_ph: raw voltage or voltage counts
    :param temperature: Temperature in degrees C
    :param salinity: Salinity in PSU
    :param pressure: Pressure in dbar
    :param coefs: External pH coefficients
    :param ph_units: The units for raw_ph, defaults to "counts"
    :param formula_version: The version of the pH formula, where
        "legacy" refers to the formula used by Fathom v3.0.4 and UCI
        v4.0.x, and "1.3" refers to the version of the Argo pH doc in
        the description
    :return: Total external pH
    """
    if ph_counts is not None:
        warnings.warn("Deprecated, use raw_ph", DeprecationWarning)
        ph_volts = convert_ph_voltage_counts(ph_counts)
    elif ph_units == 'counts':
        ph_volts = convert_ph_voltage_counts(raw_ph)
    else:  # ph_counts == 'volts'
        ph_volts = raw_ph

    t_kelvin = temperature + KELVIN_OFFSET_0C
    p_bar = pressure / 10
    f_p = _pressure_response(pressure, coefs)
    nernst = _calculate_nernst(t_kelvin)
    s_t = _total_sulfate_in_seawater(salinity)
    k_stp = _acid_dissociation_constant_of_hso4_tp(salinity, temperature, p_bar)

    term_2 = np.log10(_total_chloride_in_seawater(salinity))
    term_4 = np.log10(1 - 1.005e-3 * salinity)
    term_5 = np.log10(1 + s_t / k_stp)

    if coefs.k2_poly_order == 0:
        k2_poly = Polynomial([coefs.k2])
    else:
        k2_poly = Polynomial([coefs.k2f0, coefs.k2f1, coefs.k2f2, coefs.k2f3])
    
    eot = k2_poly(pressure) * temperature

    if formula_version == 'legacy':
        thermal_pressure = _calculate_thermal_pressure(temperature, p_bar)
        term_1 = (ph_volts - coefs.k0 - eot - f_p - thermal_pressure) / nernst
        term_3 = 2 * _log_of_hcl_activity_coefficient_of_t(salinity, temperature)

    elif formula_version == '1.3':
        term_1 = (ph_volts - coefs.k0 - eot - f_p) / nernst
        term_3 = 2 * _log_of_hcl_activity_coefficient_of_tp(salinity, temperature, p_bar)

    ph = term_1 + term_2 + term_3 - term_4 - term_5
    return ph


def convert_seafet_temperature(raw_temp, coefs: cc.TemperatureSeaFETCoefficients):
    """Converts the raw SeaFET temperature value to ITS-90 Celsius.

    :param raw_temp: raw temperature values
    :return: ITS-90 Celsius.
    """
    temp_log = np.log(raw_temp)

    temp = 1 / (
        ((coefs.tdfa3 * temp_log + coefs.tdfa2) * temp_log + coefs.tdfa1) * temp_log + coefs.tdfa0
    )

    temp_c = temp - KELVIN_OFFSET_0C

    return temp_c


def convert_internal_seafet_temperature(temperature_counts: np.ndarray):
    """Converts the raw internal temperature counts to degrees Celsius

    :param temperature_counts: raw internal temperature counts
    :return: internal temperature in Celsius
    """
    slope = 175.72
    offset = -46.85
    int_16bit = 2**16
    temperature = temperature_counts / int_16bit * slope + offset
    return temperature


def convert_seafet_relative_humidity(humidity_counts: np.ndarray, temperature: np.ndarray):
    """Convert relative humidity counts to percent

    :param humidity_counts: raw relative humidity counts
    :param temperature: converted internal temperature in Celsius
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


def convert_altimeter(
    volts: np.ndarray,
    coefs: cc.AltimeterCoefficients,
):
    """Converts a raw voltage value for altimeter.

    All equation information comes from application note 95

    :param volts: raw output voltage from altimeter sensor
    :param coefs: slope and offset for the altimeter sensors

    :return: converted height in meters
    """
    ALTIMETER_SCALAR = 300

    height = ALTIMETER_SCALAR * volts / coefs.slope - coefs.offset

    return height
