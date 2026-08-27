"""A collection of raw data conversion functions."""

import math
from collections.abc import Callable
from typing import Literal

import gsw
import numpy as np
from numpy.polynomial import Polynomial
from scipy import stats
import seawater as sw

import seabirdscientific.cal_coefficients as cc
import seabirdscientific.constants as const
import seabirdscientific.instrument_data as si
from seabirdscientific import eos80_conversion as eos80


def convert_temperature_units(
    temperature: np.ndarray,
    from_standard: Literal["ITS90", "IPTS68"],
    from_units: Literal["C", "F"],
    to_standard: Literal["ITS90", "IPTS68"] = "ITS90",
    to_units: Literal["C", "F"] = "C",
) -> np.ndarray:
    """Convert temperature between different units and calibration standards.

    Converts temperature values between Celsius and Fahrenheit, and between
    ITS90 and IPTS68 calibration standards.

    :param temperature: Temperature values to convert
    :param from_standard: Input calibration standard (ITS90 or IPTS68)
    :param from_units: Input temperature units (C or F)
    :param to_standard: Output calibration standard, defaults to ITS90
    :param to_units: Output temperature units, defaults to C

    :return: Temperature values converted to specified units and standard
    """
    _temperature = temperature.copy()
    if from_standard == "ITS90" and to_standard == "IPTS68":
        _temperature *= const.ITS90_TO_IPTS68
    elif from_standard == "IPTS68" and to_standard == "ITS90":
        _temperature /= const.ITS90_TO_IPTS68

    if from_units == "F" and to_units == "C":
        _temperature = (_temperature - 32) * 5 / 9
    if from_units == "C" and to_units == "F":
        _temperature = _temperature * 9 / 5 + 32

    return _temperature


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
    :param standard: whether to convert to ITS90 or IPTS-68 calibration
        standard
    :param units: whether to convert to celsius or fahrenheit
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
    ) - const.KELVIN_OFFSET_0C

    temperature = convert_temperature_units(temperature, "ITS90", "C", standard, units)

    return temperature


def convert_temperature_frequency(
    frequency: np.ndarray,
    coefs: cc.TemperatureFrequencyCoefficients,
    standard: Literal["ITS90", "IPTS68"] = "ITS90",
    units: Literal["C", "F"] = "C",
):
    """Convert raw frequency to temperature in degrees Celsius or
    degrees Fahrenheit

    :param frequency: raw frequency from the temperature sensor
    :param coefs: calibration coefficients for the temperature sensor
    :return: temperature in Celsius or Fahrenheit
    """
    fLog = np.log(coefs.f0 / frequency)
    temperature = (
        1 / (coefs.g + coefs.h * fLog + coefs.i * fLog**2 + coefs.j * fLog**3)
        - const.KELVIN_OFFSET_0C
    )

    temperature = convert_temperature_units(temperature, "ITS90", "C", standard, units)

    return temperature


def convert_pressure_units(
    pressure: np.ndarray,
    from_units: Literal["dbar", "psia", "psig"],
    to_units: Literal["dbar", "psia", "psig"] = "psia",
) -> np.ndarray:
    """Convert pressure between different units.

    Converts pressure values between decibars (dbar), pounds per square inch
    absolute (psia), and pounds per square inch gauge (psig).

    :param pressure: Pressure values to convert
    :param from_units: Input pressure units (dbar, psia, or psig)
    :param to_units: Output pressure units, defaults to psia

    :return: Pressure values converted to specified units
    """
    _pressure = pressure.copy()

    if from_units == "psia" and to_units in ("dbar", "psig"):
        _pressure -= const.SEA_LEVEL_PRESSURE
    elif from_units in ("dbar", "psig") and to_units == "psia":
        _pressure += const.SEA_LEVEL_PRESSURE

    if from_units in ("psia", "psig") and to_units == "dbar":
        _pressure *= const.PSI_TO_DBAR
    elif from_units == "dbar" and to_units in ("psia", "psig"):
        _pressure /= const.PSI_TO_DBAR

    return _pressure


def convert_pressure(
    pressure_count: np.ndarray,
    compensation_voltage: np.ndarray,
    coefs: cc.PressureCoefficients,
    units: Literal["dbar", "psia", "psig"] = "psia",
):
    """Converts pressure counts to sea pressure (psig and dbar) and
    absolute pressure (psia)

    pressure_count and compensation_voltage are expected to be raw data
    from an instrument in A/D counts

    :param pressure_count: pressure value to convert, in A/D counts
    :param compensation_voltage: pressure temperature compensation
        voltage, in counts or volts depending on the instrument
    :param coefs: calibration coefficients for the pressure sensor
    :param units: whether or not to use dbar, psig, or psia as the
        returned unit type, defaults to psia

    :return: sea pressure val in dbar, psig, or psia according to units
    """

    t = (
        coefs.ptempa0
        + coefs.ptempa1 * compensation_voltage
        + coefs.ptempa2 * compensation_voltage**2
    )
    x = pressure_count - coefs.ptca0 - coefs.ptca1 * t - coefs.ptca2 * t**2
    n = x * coefs.ptcb0 / (coefs.ptcb0 + coefs.ptcb1 * t + coefs.ptcb2 * t**2)
    pressure = coefs.pa0 + coefs.pa1 * n + coefs.pa2 * n**2

    pressure = convert_pressure_units(pressure, "psia", units)

    return pressure


def convert_pressure_digiquartz(
    pressure_count: np.ndarray,
    compensation_voltage: np.ndarray,
    coefs: cc.PressureDigiquartzCoefficients,
    units: Literal["dbar", "psia", "psig"],
    sample_interval: float,
):
    """Converts pressure counts to PSIA (pounds per square inch,
    abolute), PSIG, or dbar for a digiquartz pressure sensor.

    :param pressure_count: pressure value to convert, in A/D counts
    :param compensation_voltage: pressure temperature, in A/D counts
    :param coefs: calibration coefficients for the digiquartz pressure
        sensor
    :param units: whether or not to use dbar, psig, or psia as the
        returned unit type, defaults to psia
    :param sample_interval: sample rate of the data to be used for
        temperature compensation correction, in seconds
    :return: pressure val in dbar, psig, or psia according to units
    """

    # First, average temperature compensation over 30 seconds
    def modification_function(x):
        return x * coefs.ad590m + coefs.ad590b

    # using a short name to make the equations a little easier to read
    v = _compute_rolling_average(compensation_voltage, 30, sample_interval, modification_function)

    # Now, calculate pressure
    t = 1 / pressure_count * 1e6  # convert to period in usec
    c = coefs.c1 + coefs.c2 * v + coefs.c3 * v**2
    d = coefs.d1 + coefs.d2 * v
    t0 = coefs.t1 + coefs.t2 * v + coefs.t3 * v**2 + coefs.t4 * v**3 + coefs.t5 * v**4

    one_minus_t_ratio = 1 - (t0**2) / (t**2)
    # p is absolute pressure according to Paroscientific cal sheet
    p = c * one_minus_t_ratio * (1 - d * one_minus_t_ratio)

    pressure = convert_pressure_units(p, "psia", units)

    return pressure


def convert_conductivity_units(
    conductivity: np.ndarray,
    from_units: Literal["S/m", "mS/cm", "uS/cm"],
    to_units: Literal["S/m", "mS/cm", "uS/cm"] = "S/m",
) -> np.ndarray:
    """Convert conductivity between different units.

    Converts conductivity values between Siemens per meter (S/m),
    milliSiemens per centimeter (mS/cm), and microSiemens per centimeter (uS/cm).

    :param conductivity: Conductivity values to convert
    :param from_units: Input conductivity units (S/m, mS/cm, or uS/cm)
    :param to_units: Output conductivity units, defaults to S/m

    :return: Conductivity values converted to specified units
    """
    if from_units == "S/m" and to_units == "mS/cm":
        return conductivity * 10
    elif from_units == "S/m" and to_units == "uS/cm":
        return conductivity * 1e4
    elif from_units == "mS/cm" and to_units == "S/m":
        return conductivity / 10
    elif from_units == "mS/cm" and to_units == "uS/cm":
        return conductivity * 1000
    elif from_units == "uS/cm" and to_units == "S/m":
        return conductivity / 1e4
    elif from_units == "uS/cm" and to_units == "mS/cm":
        return conductivity / 1000
    else:
        return conductivity


def convert_conductivity(
    conductivity_count: np.ndarray,
    temperature: np.ndarray,
    pressure: np.ndarray,
    coefs: cc.ConductivityCoefficients,
    instrument_type: si.InstrumentType,
    units: Literal["S/m", "mS/cm", "uS/cm"] = "S/m",
):
    """Converts raw conductivity counts to S/m, mS/cm, or uS/cm.

    Data is expected to be raw data from instrument in A/D counts

    :param conductivity_count: conductivity value to convert, in A/D
        counts
    :param temperature: reference temperature, in degrees C
    :param pressure: reference pressure, in dbar
    :param coefs: calibration coefficient for the conductivity sensor
    :param instrument_type: the instrument that recorded conductivity
    :param units: the conductivity units to convert to, defaults to S/m

    :return: conductivity val converted to S/m
    """
    scalar = 1
    if instrument_type == si.InstrumentType.SBE911Plus:
        scalar = 1 / 10

    if instrument_type in (
        si.InstrumentType.SBE16Plus,
        si.InstrumentType.SBE19Plus,
        si.InstrumentType.SBE911Plus,
    ):
        f = conductivity_count / 1000
    else:
        f = conductivity_count * np.sqrt(1 + coefs.wbotc * temperature) / 1000

    numerator = coefs.g + coefs.h * f**2 + coefs.i * f**3 + coefs.j * f**4
    denominator = 1 + coefs.ctcor * temperature + coefs.cpcor * pressure
    conductivity = convert_conductivity_units(numerator / denominator * scalar, "S/m", units)
    return conductivity


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
    water_type: Literal["salt", "fresh"] = "salt",
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
        pressure /= const.DBAR_TO_PSI

    if water_type == "fresh":
        depth = pressure * const.FRESHWATER_PRESSURE_TO_DEPTH
    else:
        depth = -gsw.z_from_p(pressure, latitude)

    if depth_units == "ft":
        depth *= const.METERS_TO_FEET

    return depth




def convert_sbe63_oxygen(
    raw_oxygen_phase: np.ndarray,
    thermistor: np.ndarray,
    pressure: np.ndarray,
    salinity: np.ndarray,
    coefs: cc.Oxygen63Coefficients,
    thermistor_coefs: cc.Thermistor63Coefficients,
    thermistor_units: Literal["volts", "C"] = "volts",  # Is this volts or frequency?
    units: Literal[
        "ml/l",
        "mg/l",
        "umol/kg",
        "umol/l",
        "saturation_percent",
        "ox_temperature_c",
        "ox_temperature_f",
        "raw_phase_usec",
        "raw_phase_v",
    ] = "ml/l",
    external_temperature: np.ndarray | None = None,
):
    """Returns the data after converting it to desired units.

    raw_oxygen_phase is expected to be in raw phase, raw_thermistor_temp
    in counts, pressure in dbar, and salinity in practical salinity (PSU)

    :param raw_oxygen_phase: SBE63 phase value, in microseconds
    :param thermistor_temp: SBE63 thermistor data to use are reference,
        in volts or degrees C (see thermistor_units param)
    :param pressure: Converted pressure value from the attached CTD, in
        dbar
    :param salinity: Converted salinity value from the attached CTD, in
        practical salinity PSU
    :param coefs (cc.Oxygen63Coefficients): calibration coefficients for
        the SBE63 sensor
    :param thermistor_coefs (cc.Thermistor63Coefficients): calibration coefficients for
        the SBE63 thermistor sensor
    :param thermistor_units: units of thermistor_temp input
    :param units: the units to return the oxygen values in. Options are:
        ml/l, mg/l, umol/kg, umol/l, saturation_percent, ox_temperature_c, ox_temperature_f, raw_phase_usec, raw_phase_v. Defaults to ml/l.
    :param external_temperature: optional external temperature to use for oxygen conversion, in degrees C. Required for umol/kg and percentage saturation units. If not provided, the thermistor will be used for temperature.

    :return: converted Oxygen value, in ml/l
    """
    if thermistor_units == "volts":
        temperature = convert_sbe63_thermistor(thermistor, thermistor_coefs)
    elif thermistor_units == "C":
        temperature = thermistor
    else:
        raise ValueError

    if units == "ox_temperature_c":
        return temperature
    elif units == "ox_temperature_f":
        return temperature * 9 / 5 + 32  # Convert C to F
    elif units == "raw_phase_usec":
        return raw_oxygen_phase

    oxygen_volts = raw_oxygen_phase / const.OXYGEN_PHASE_TO_VOLTS  # from the manual

    if units == "raw_phase_v":
        return oxygen_volts

    ksv = coefs.c0 + coefs.c1 * temperature + coefs.c2 * temperature**2

    s_corr_exp = _compute_ln_salinity_correction(temperature, salinity)
    s_corr = math.e**s_corr_exp

    # temperature in Kelvin
    temperature_k = temperature + const.KELVIN_OFFSET_0C
    p_corr_exp = (coefs.e * pressure) / temperature_k
    p_corr = math.e**p_corr_exp

    # fmt: off
    oxygen = (
        (((coefs.a0 + coefs.a1 * temperature + coefs.a2 * oxygen_volts**2)
        / (coefs.b0 + coefs.b1 * oxygen_volts) - 1.0) / ksv) * s_corr * p_corr
    )
    # fmt: on
    # If an external temperature is provided, use that for the gsw functions instead of thermistor
    temperature_to_use = external_temperature if external_temperature is not None else temperature
    if units == "ml/l":
        return oxygen
    elif units == "mg/l":
        return convert_oxygen_to_mg_per_l(oxygen)
    elif units == "umol/kg":
        potential_density = potential_density_from_t_s_p(temperature_to_use, salinity, pressure)
        return convert_oxygen_to_umol_per_kg(oxygen, potential_density)
    elif units == "umol/l":
        return convert_oxygen_to_umol_per_l(oxygen)
    elif units == "saturation_percent":
        # O2 Saturation always uses GG calc, as it is more accurate than Weiss
        oxygen_saturation = derive_oxygen_saturation_gg(temperature_to_use, salinity)
        oxygen_saturation_percent = oxygen * 100 / oxygen_saturation

        # handle cases where oxygen saturation is flagged
        return np.where(
            oxygen_saturation != const.FLAG_VALUE, oxygen_saturation_percent, const.FLAG_VALUE
        )


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
        - const.KELVIN_OFFSET_0C
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
    units: Literal[
        "ml/l", "mg/l", "umol/kg", "umol/l", "dov/dt", "saturation_percent", "raw_voltage"
    ] = "ml/l",
):
    """Returns the data after converting it to desired units.

    voltage is expected to be in volts, temperature in ITS-90 deg c, pressure
    in dbar, and salinity in practical salinity (PSU). All equation
    information comes from Application Note 64

    :param voltage: SBE43 voltage
    :param temperature: temperature value converted to ITS-90 deg C
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
    :param units: the units to return the oxygen values in. Options are:
        ml/l, mg/l, umol/kg, umol/l, dov/dt, saturation_percent, raw_voltage
        Defaults to ml/l.

    :return: converted Oxygen values, in ml/l
    """
    # start with all 0 for the dvdt
    dvdt_values = np.zeros(len(voltage))
    if apply_tau_correction or units == "dov/dt":
        # Calculates how many scans to have on either side of our median
        # point, accounting for going out of index bounds
        scans_per_side = math.floor(window_size / 2 / sample_interval)
        for i in range(scans_per_side, len(voltage) - scans_per_side):
            ox_subset = voltage[i - scans_per_side : i + scans_per_side + 1]

            time_subset = np.arange(
                0, len(ox_subset) * sample_interval, sample_interval, dtype=float
            )

            result = stats.linregress(time_subset, ox_subset)

            dvdt_values[i] = result.slope

    if units == "dov/dt":
        return dvdt_values

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

    if units == "raw_voltage":
        # Return the corrected voltage values if the user wants raw voltage
        return correct_ox_voltages

    oxygen = _convert_sbe43_oxygen(
        correct_ox_voltages,
        temperature,
        pressure,
        salinity,
        coefs,
        dvdt_values,
    )
    if units == "ml/l":
        return oxygen
    elif units == "mg/l":
        return convert_oxygen_to_mg_per_l(oxygen)
    elif units == "umol/kg":
        potential_density = potential_density_from_t_s_p(temperature, salinity, pressure)
        return convert_oxygen_to_umol_per_kg(oxygen, potential_density)
    elif units == "umol/l":
        return convert_oxygen_to_umol_per_l(oxygen)
    elif units == "saturation_percent":
        # O2 Saturation always uses GG calc, as it is more accurate than Weiss
        oxygen_saturation = derive_oxygen_saturation_gg(temperature, salinity)
        oxygen_saturation_percent = oxygen * 100 / oxygen_saturation
        # handle cases where oxygen saturation is flagged
        return np.where(
            oxygen_saturation != const.FLAG_VALUE, oxygen_saturation_percent, const.FLAG_VALUE
        )


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
    information comes from Application Note 64.
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

    # Oxygen Solubility equation constants, From Application Note 64 Appendix A
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

    ts = _compute_scaled_temperature(temperature)
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
        * np.exp((coefs.e * pressure) / (temperature + const.KELVIN_OFFSET_0C))
    )
    return oxygen


def convert_oxygen_units(
    oxygen: np.ndarray,
    temperature: np.ndarray,
    pressure: np.ndarray,
    salinity: np.ndarray,
    from_units: Literal["ml/l", "mg/l", "umol/kg", "umol/l", "saturation_percent"],
    to_units: Literal["ml/l", "mg/l", "umol/kg", "umol/l", "saturation_percent"],
):
    """Convert oxygen values between supported oxygen units.

    Conversion is always done in two steps:
    1) convert input values to ml/L
    2) convert ml/L values to the target units

    :param oxygen: oxygen values in ``from_units``
    :param temperature: temperature in degrees C (used for saturation and density)
    :param pressure: pressure in dbar (used for density)
    :param salinity: salinity in PSU (used for saturation and density)
    :param from_units: source oxygen units
    :param to_units: destination oxygen units

    :return: oxygen values in ``to_units``
    """
    if from_units == to_units:
        return oxygen.copy()

    potential_density = None

    # Step 1: normalize to ml/l
    if from_units == "ml/l":
        oxygen_ml_per_l = oxygen
    elif from_units == "mg/l":
        oxygen_ml_per_l = oxygen / const.OXYGEN_MLPERL_TO_MGPERL
    elif from_units == "umol/l":
        oxygen_ml_per_l = oxygen / const.OXYGEN_MLPERL_TO_UMOLPERL
    elif from_units == "umol/kg":
        potential_density = potential_density_from_t_s_p(temperature, salinity, pressure)
        oxygen_ml_per_l = oxygen * (potential_density + 1000) / const.OXYGEN_MLPERL_TO_UMOLPERKG
    elif from_units == "saturation_percent":
        oxygen_saturation = derive_oxygen_saturation_gg(temperature, salinity)
        oxygen_ml_per_l = oxygen_saturation * oxygen / 100.0
    else:
        raise ValueError(f"Unsupported from_units: {from_units}")

    # Step 2: convert from ml/l to target units
    if to_units == "ml/l":
        converted = oxygen_ml_per_l
    elif to_units == "mg/l":
        converted = convert_oxygen_to_mg_per_l(oxygen_ml_per_l)
    elif to_units == "umol/l":
        converted = convert_oxygen_to_umol_per_l(oxygen_ml_per_l)
    elif to_units == "umol/kg":
        if potential_density is None:
            potential_density = potential_density_from_t_s_p(temperature, salinity, pressure)
        converted = convert_oxygen_to_umol_per_kg(oxygen_ml_per_l, potential_density)
    elif to_units == "saturation_percent":
        oxygen_saturation = derive_oxygen_saturation_gg(temperature, salinity)
        oxygen_saturation_percent = oxygen_ml_per_l * 100 / oxygen_saturation
        converted = np.where(
            oxygen_saturation != const.FLAG_VALUE,
            oxygen_saturation_percent,
            const.FLAG_VALUE,
        )
    else:
        raise ValueError(f"Unsupported to_units: {to_units}")

    # Preserve explicit input bad flags where present.
    return np.where(oxygen == const.FLAG_VALUE, const.FLAG_VALUE, converted)


def convert_oxygen_to_mg_per_l(ox_values: np.ndarray):
    """Converts given oxygen values to milligrams/Liter.

    From Application Note 64.

    :param ox_values: oxygen values, already converted to ml/L

    :return: oxygen values converted to milligrams/Liter
    """

    return ox_values * const.OXYGEN_MLPERL_TO_MGPERL


def convert_oxygen_to_umol_per_kg(ox_values: np.ndarray, potential_density: np.ndarray):
    """Converts given oxygen values to micromoles/kg.

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

    :return: oxygen values converted to micromoles/kg
    """

    oxygen_umolkg = (ox_values * const.OXYGEN_MLPERL_TO_UMOLPERKG) / (potential_density + 1000)
    return oxygen_umolkg


def convert_oxygen_to_umol_per_l(ox_values: np.ndarray):
    """Converts given oxygen values to micromoles/l.
    :param ox_values: oxygen values, already converted to ml/L

    :return: oxygen values converted to micromoles/ml
    """

    oxygen_umolkg = ox_values * const.OXYGEN_MLPERL_TO_UMOLPERL
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
        1.98416e-4 * (temperature + const.KELVIN_OFFSET_0C) * coefs.slope
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
        nitrate *= const.UMNO3_TO_MGNL

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
    nernst_term = const.R * temperature * np.log(10) / const.F
    return nernst_term


def convert_internal_seafet_ph(
    raw_ph: np.ndarray,
    temperature: np.ndarray,
    coefs: cc.PHSeaFETInternalCoefficients,
    ph_units: Literal["counts", "volts"] = "counts",
):
    """Calculates the internal pH on the total scale given the
    temperature and internal FET voltage

    :param raw_ph: Raw voltage or voltage counts
    :param temperature: Sample temperature
    :param coefs: SeaFET calibration coefficients
    :param ph_units: The units of raw_ph, defaults to 'counts'
    :return: calculated pH on the total scale for the SeaFET internal
        reference
    """
    if ph_units == "counts":
        ph_volts = convert_ph_voltage_counts(raw_ph)
    else:  # ph_counts == 'volts'
        ph_volts = raw_ph

    nernst = _calculate_nernst(temperature + const.KELVIN_OFFSET_0C)
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
    t_kelvin = temperature + const.KELVIN_OFFSET_0C
    term_2 = (v_hcl * pressure) / (np.log(10) * const.R * t_kelvin * 10) / 2
    log_y_hcl_tp = log_y_hcl + term_2
    return log_y_hcl_tp


def _acid_dissociation_constant_of_hso4(salinity: np.ndarray, temperature: np.ndarray):
    """From SBS application note 99. Calculated as (Dickson et al. 2007)

    :param salinity: Salinity in PSU
    :param temperature: Temperature in Kelvin
    :return: Acid dissociation constant of HSO4
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
    t_kelvin = temperature + const.KELVIN_OFFSET_0C
    k_s = _acid_dissociation_constant_of_hso4(salinity, t_kelvin)
    v_bar_s = _partial_molal_hso4_volume(temperature)
    k_bar_s = _hso4_compressibility(temperature)
    exponent = (-v_bar_s * pressure + 0.5 * k_bar_s * pressure**2) / (const.R * t_kelvin * 10)
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
    raw_ph: np.ndarray,
    temperature: np.ndarray,
    salinity: np.ndarray,
    pressure: np.ndarray,
    coefs: cc.PHSeaFETExternalCoefficients,
    ph_units: Literal["counts", "volts"] = "counts",
    formula_version: Literal["legacy", "1.3"] = "1.3",
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
    if ph_units == "counts":
        ph_volts = convert_ph_voltage_counts(raw_ph)
    else:  # ph_counts == 'volts'
        ph_volts = raw_ph

    t_kelvin = temperature + const.KELVIN_OFFSET_0C
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

    if formula_version == "legacy":
        thermal_pressure = _calculate_thermal_pressure(temperature, p_bar)
        term_1 = (ph_volts - coefs.k0 - eot - f_p - thermal_pressure) / nernst
        term_3 = 2 * _log_of_hcl_activity_coefficient_of_t(salinity, temperature)

    elif formula_version == "1.3":
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

    temp_c = temp - const.KELVIN_OFFSET_0C

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


def buoyancy_frequency(
    temperature: np.ndarray,
    salinity: np.ndarray,
    pressure: np.ndarray,
    gravity: float,
):
    """Calculates an N^2 value (buoyancy frequency) for the given window
    of temperature, salinity, and pressure, at the given latitude.

    Expect temperature as conservative temperature, salinity as abslute
    salinity, and pressure as dbar, all of the same length. Performs the
    calculation using TEOS-10 and specific volume.

    :param temperature: temperature values for the given window
    :param salinity: salinity values for the given window
    :param pressure: pressure values for the given window
    :param gravity: gravity value

    :return: A single N^2 [Brunt-Väisälä (buoyancy) frequency]
    """

    db_to_pa = 1e4
    # Wrap these as a length-1 array so that GSW accepts them
    mean_pressure = [np.mean(pressure)]
    mean_temperature = [np.mean(temperature)]
    mean_salinity = [np.mean(salinity)]

    # Compute average specific volume, temp expansion ceoff,
    # and saline contraction coeff over window
    (specific_volume, alpha, beta) = gsw.specvol_alpha_beta(
        mean_salinity, mean_temperature, mean_pressure
    )

    # Estimate vertical gradient of conservative temp
    dct_dp = stats.linregress(pressure, temperature)
    # TODO: error handling with r, p, std_error

    # Estimate vertical gradient of absolute salinity
    dsa_dp = stats.linregress(pressure, salinity)
    # TODO: error handling with r, p, std_error

    # Compute N2 combining computed ceofficients and vertical gradients.
    # we index into specific_volume, alpha, and beta as they are all arrays of len 1
    n2 = gravity**2 / (specific_volume[0] * db_to_pa)
    n2 *= beta[0] * dsa_dp.slope - alpha[0] * dct_dp.slope
    return n2


def buoyancy(
    temperature: np.ndarray,
    salinity: np.ndarray,
    pressure: np.ndarray,
    latitude: np.ndarray,
    longitude: np.ndarray,
    window_size: float,
    use_modern_formula=True,
    flag_value=const.FLAG_VALUE,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Calculates the 4 buoyancy values based off the incoming data.

    Data is expected to have already been binned via Bin_Average using
    decibar pressure bins. All arrays are expected to be the same
    length, except for latitude and longitude, which can be length 1.
    Optionally can use the former calculation for buoyancy frequency
    from the SBE Data Processing Manual, but defaults to a newer formula
    using TEOS-10.

    :param temperature_c: Temperature in ITS-90 degrees C
    :param salinity_prac: Practical salinity in PSU
    :param pressure_dbar: Pressure in dbar
    :param latitude: latitude values. If length 1, gets applied to all
        values.
    :param longitude: longitude values. If length 1, gets applied to all
        values.
    :param window_size: window size to use. If this number is smaller
        than the binned window size, round up to a minium of 3 scans.
        I.E. uses the center scan and one scan on each side of it at the
        very least
    :param use_modern_formula: Whether to use a modern formula for
        calculating buoyancy frequency. Defaults to true.
    :param flag_value: Bad Flag value to use for marking bad scans.
        Defaults to -9.99e-29

    :return: a tuple of ndarrays including: buoyancy frequency squared,
        buoyancy frequency, stability, and scaled stability
    """

    _salinity, _temperature, _pressure, _latitude, _longitude = np.broadcast_arrays(
        salinity, temperature, pressure, latitude, longitude
    )

    # Get the original bin size that we're working with, using the
    # second and third bin so we don't have to worry about the surface
    # bin
    original_bin_size = abs(_pressure[2] - _pressure[1])

    # Calculates how many scans to have on either side of our median
    # point, but need at least 1 (for a total of 3 scans)
    scans_per_side = max(math.floor(window_size / original_bin_size / 2), 1)

    salinity_abs = gsw.SA_from_SP(_salinity, _pressure, _longitude, _latitude)
    temperature_conservative = gsw.CT_from_t(salinity_abs, _temperature, _pressure)

    # create our result np.ndarrays with the flag value as default
    buoyancy_freq_squared = np.full(len(_temperature), flag_value)
    buoyancy_freq = np.full(len(_temperature), flag_value)
    stability = np.full(len(_temperature), flag_value)
    scaled_stability = np.full(len(_temperature), flag_value)

    # start loop at scans_per_side
    for i in range(scans_per_side, len(temperature_conservative) - scans_per_side):
        min_index = i - scans_per_side
        max_index = (
            i + scans_per_side + 1
        )  # add + 1 because slicing does not include the max_index

        pressure_subset = _pressure[min_index:max_index]
        temperature_cons_subset = temperature_conservative[min_index:max_index]
        temperature_its_subset = _temperature[min_index:max_index]
        salinity_subset = salinity_abs[min_index:max_index]

        mean_pressure = [np.mean(pressure_subset)]
        gravity = gsw.grav([_latitude[i]], mean_pressure)[0]

        if use_modern_formula:
            salinity_subset = salinity_abs[min_index:max_index]
            n2 = buoyancy_frequency(
                temperature_cons_subset, salinity_subset, pressure_subset, gravity
            )
        else:
            salinity_subset = _salinity[min_index:max_index]
            n2 = eos80.bouyancy_frequency(
                temperature_its_subset, salinity_subset, pressure_subset, gravity
            )

        buoyancy_freq_squared[i] = n2
        if n2 >= 0:
            buoyancy_freq[i] = math.sqrt(n2) * 3600 / (2 * np.pi)
        else:
            # using the negative square root of the absolute buoyancy squared to match seasoft
            buoyancy_freq[i] = -math.sqrt(abs(n2)) * 3600 / (2 * np.pi)
        stability[i] = n2 / gravity
        scaled_stability[i] = stability[i] * 1e8

    return (buoyancy_freq_squared, buoyancy_freq, stability, scaled_stability)


def _compute_scaled_temperature(temperature: np.ndarray) -> np.ndarray:
    return np.log((const.KELVIN_OFFSET_25C - temperature) / (const.KELVIN_OFFSET_0C + temperature))


def _compute_ln_salinity_correction(temperature: np.ndarray, salinity: np.ndarray) -> np.ndarray:
    """Compute natural logarithm of the salinity correction for Garcia and Gordon
    Oxygen Solubility. Also applicable to SBE 63 Oxygen

    :param temperature: Temperature in degrees Celsius
    :param salinity: Salinity in PSU

    :return: Natural logarithm of the salinity correction"""
    sol_b0 = -6.24523e-3
    sol_b1 = -7.37614e-3
    sol_b2 = -1.0341e-2
    sol_b3 = -8.17083e-3
    sol_c0 = -4.88682e-7

    ts = _compute_scaled_temperature(temperature)
    s_corr = (
        salinity * (sol_b0 + sol_b1 * ts + sol_b2 * ts**2 + sol_b3 * ts**3) + sol_c0 * salinity**2
    )
    return s_corr


def _compute_rolling_average(
    compute_var: np.ndarray,
    window_size: float,
    sample_interval: float,
    modification_fn: Callable | None = None,
) -> np.ndarray:
    """Computes a rolling average of the given variable over the specified window size.
    Averages with equal number values on either side of the center of each window.

    :param compute_var: The variable to compute the rolling average for.
    :param window_size: The size of the rolling window (in seconds).
    :param sample_interval: The time interval between samples (in seconds).
    :param modification_fn: Optional function to modify the computed rolling average.

    :return: An array containing the rolling average values.
    """
    if window_size <= 0:
        raise ValueError("Window size must be a positive integer.")

    # Calculate the number of samples in the rolling window
    num_samples = int(window_size / sample_interval)

    # Determine padding needed for both ends
    pad_before = num_samples // 2
    pad_after = num_samples - 1 - pad_before

    # Pad the array using the edge values
    # This prevents the ends from dropping off or pulling toward zero
    padded_data = np.pad(compute_var, (pad_before, pad_after), mode="edge")

    # Compute the rolling average using numpy's convolve function
    weights = np.ones(num_samples) / num_samples
    rolling_avg = np.convolve(padded_data, weights, mode="valid")

    # Apply the modification function if provided
    if modification_fn is not None:
        rolling_avg = modification_fn(rolling_avg)

    return rolling_avg


def derive_descent_rate(
    depth: np.ndarray,
    window_size: float,
    sample_interval: float,
) -> np.ndarray:
    """Derives the descent rate from the depth values.

    :param depth: Depth values in meters or feet.
    :param window_size: Window size to use for the derivative calculation in seconds
    :param sample_interval: Sample interval in seconds

    :return: np.ndarray of descent rate values in meters per second or feet per second, depending on the input depth units.
    """
    # TODO: slightly different calculation from sbe data processing, but this is was more simple

    # Calculate the number of samples to include in the window based on the sample interval
    samples_per_window = max(int(window_size / sample_interval + 1), 1)
    samples_per_side = max(int(samples_per_window // 2), 1)
    time_array = np.arange(len(depth)) * sample_interval

    # Calculate the descent rate using a centered difference method
    descent_rate = np.full(len(depth), 0.0)  # Initialize with 0 for edge cases

    for i in range(samples_per_side, len(depth) - samples_per_side):
        # linear regression for descent rate on subset of depth and time
        slope, _, _, _, _ = stats.linregress(
            time_array[i - samples_per_side : i + samples_per_side + 1],
            depth[i - samples_per_side : i + samples_per_side + 1],
        )
        descent_rate[i] = slope

    return descent_rate


def derive_acceleration(
    depth: np.ndarray,
    window_size: float,
    sample_interval: float,
) -> np.ndarray:
    """Derives the acceleration from the depth values.

    :param depth: Depth values in meters or feet.
    :param window_size: Window size to use for the derivative calculation in seconds
    :param sample_interval: Sample interval in seconds

    :return: np.ndarray of acceleration values in meters per second squared or feet per second squared, depending on the input depth units.
    """
    # Calculate the number of samples to include in the window based on the sample interval
    descent_rate = derive_descent_rate(depth, window_size, sample_interval)

    # Calculate the acceleration using a centered difference method
    acceleration = np.full(len(depth), 0.0)  # Initialize with 0 for edge cases

    for i in range(1, len(depth)):
        # Follow SBE Processing calc: Acc = (DescentRate[i] - DescentRate[i-1]) / SampleInterval
        acceleration[i] = (descent_rate[i] - descent_rate[i - 1]) / sample_interval

    return acceleration


def derive_oxygen_saturation_gg(
    temperature: np.ndarray,
    salinity: np.ndarray,
    flag_value=const.FLAG_VALUE,
):
    """Calculates the oxygen saturation in ml/L.

    From Garcia and Gordon L&O 37(6), 1992, 1307 - 1312
    Provide better fit and better estimation of o2 solubility at end members
    Note: SBE Data Processing returns -99 for t < -5, t > 50, s < 0, and s > 60
    This software sets these to flag value instead of -99

    :param temperature: temperature in degrees C
    :param salinity: salinity in PSU

    :return: oxygen saturation in ml/L
    """
    ts = _compute_scaled_temperature(temperature)
    ts2 = ts**2
    ts3 = ts**3

    OA0 = 2.00907
    OA1 = 3.22014
    OA2 = 4.0501
    OA3 = 4.94457
    OA4 = -0.256847
    OA5 = 3.88767

    ox_sol = OA0 + OA1 * ts + OA2 * ts2 + OA3 * ts3 + OA4 * ts2 * ts2 + OA5 * ts2 * ts3

    ox_sol += _compute_ln_salinity_correction(temperature, salinity)
    ox_sol = np.exp(ox_sol)

    # clean up values from invalid inputs, as SBE Data Processing does
    ox_sol = np.where(temperature < -5, flag_value, ox_sol)
    ox_sol = np.where(temperature > 50, flag_value, ox_sol)
    ox_sol = np.where(salinity < 0, flag_value, ox_sol)
    ox_sol = np.where(salinity > 60, flag_value, ox_sol)

    return ox_sol


def derive_oxygen_saturation_w(
    temperature: np.ndarray,
    salinity: np.ndarray,
    flag_value=const.FLAG_VALUE,
):
    """Calculates the oxygen saturation in ml/L.

    Uses Weiss formula from 1970
    Note: SBE Data Processing returns -99 for t < 0
    This software sets these to flag value instead of -99

    :param temperature: temperature in degrees C
    :param salinity: salinity in PSU

    :return: oxygen saturation in ml/L
    """

    t0 = temperature + const.KELVIN_OFFSET_0C

    t1 = np.where(t0 > 0, 100.0 / t0, 0)
    t2 = t0 / 100.0

    A1 = -173.4292
    A2 = 249.6339
    A3 = 143.3483
    A4 = -21.8492
    B1 = -0.033096
    B2 = 0.014259
    B3 = -0.00170

    ox_sol = A1 + A2 * t1 + A3 * np.log(t2) + A4 * t2 + salinity * (B1 + B2 * t2 + B3 * t2 * t2)
    ox_sol = np.exp(ox_sol)

    # clean up values from invalid inputs, as SBE Data Processing does
    ox_sol = np.where(temperature < 0, flag_value, ox_sol)
    return ox_sol


def derive_nitrogen_saturation(
    temperature: np.ndarray,
    salinity: np.ndarray,
) -> np.ndarray:
    """Calculate nitrogen saturation from temperature and salinity.

    This is a vectorized implementation of ``N2_SaturationCalc``.

    :param temperature: temperature in degrees C
    :param salinity: salinity in PSU

    :return: nitrogen saturation values
    """

    temperature_arr, salinity_arr = np.broadcast_arrays(temperature, salinity)

    t = temperature_arr.astype(float, copy=False)
    s = salinity_arr.astype(float, copy=False)

    t0 = t + const.KELVIN_OFFSET_0C
    t1 = np.divide(100.0, t0, out=np.zeros_like(t0), where=t0 != 0.0)
    t2 = np.maximum(t0 / 100.0, 1.0e-6)

    a1 = -172.4965
    a2 = 248.4262
    a3 = 143.0738
    a4 = -21.7120
    b1 = -0.049781
    b2 = 0.025018
    b3 = -0.0034861

    n2 = a1 + a2 * t1 + a3 * np.log(t2) + a4 * t2 + s * (b1 + b2 * t2 + b3 * t2 * t2)
    n2 = np.exp(n2)

    return np.where(t0 < 0.0, 99.0, n2)


def convert_cstar_attenuation(raw: np.ndarray, coefs: cc.CstarCoefficients):
    """Converts C-Star raw voltage to attenuation.

    All equation information comes from SBE Data Processing code.

    See applicaiton note 91 for calculating m and b coefficients.

    :param raw: raw output voltage from C-Star sensor
    :param coefs: m, b, and path_length [m] for C-Star sensor

    :return: converted attenuation in %
    """
    im = (raw * coefs.m + coefs.b) / 100
    im = np.maximum(im, 0.000001)
    attenuation = (-1 / coefs.path_length) * np.log(im)

    return attenuation


def convert_cstar_transmittance(raw: np.ndarray, coefs: cc.CstarCoefficients):
    """Converts C-Star raw voltage to transmittance.

    All equation information comes from SBE Data Processing code.

    See applicaiton note 91 for calculating m and b coefficients.

    :param raw: raw output voltage from C-Star sensor
    :param coefs: m and b coefficients for C-Star sensor

    :return: converted transmittance in 1/m
    """
    transmittance = raw * coefs.m + coefs.b

    return transmittance





