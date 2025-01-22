"""A collection of classes and functions related to the processing of
instrument data.
"""

# Classes:
#   InstrumentType (Enum)
#   HexDataTypes (Enum)
#   Sensors (Enum)
#   MeasurementSeries
#   InstrumentData
# Functions:
#   cnv_to_instrument_data (Path) -> InstrumentData
#   fix_exponents (List[str]) -> List[str]
#   read_hex_file (str, InstrumentType, List[Sensors], bool) -> pd.DataFrame
#   preallocate_dataframe (InstrumentType, str, List[Sensors], bool, int) -> pd.DataFrame
#   read_hex (InstrumentType, str, List[Sensors], bool) -> dict
#   read_SBE19plus_format_0 (str, List[Sensors], bool) -> dict
#   read_SBE37SM_format_0 (str, List[Sensors]) -> dict


# Native imports
from enum import Enum
from dataclasses import dataclass
from datetime import date, datetime
from logging import getLogger
from typing import List, Dict
from pathlib import Path

# Third-party imports
import numpy as np
import pandas as pd

# Sea-Bird imports

# Internal imports

logger = getLogger(__name__)

COUNTS_TO_VOLTS = 13107
SECONDS_BETWEEN_EPOCH_AND_2000 = 946684800


class InstrumentType(Enum):
    """The type of instrument that generated the hex file being read"""

    SBE37SM = "37-SM"
    SBE37SMP = "37-SMP"
    SBE37SMPODO = "37-SMP-ODO"
    SBE37IM = "37-IM"
    SBE37IMP = "37-IMP"
    SBE37IMPODO = "37-IMP-ODO"
    SBE19Plus = "19plus" # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75


class HexDataTypes(Enum):
    """Possible data types in hex files"""

    temperature = "temperature" # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    conductivity = "conductivity" # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    pressure = "pressure" # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    temperatureCompensation = "temperature compensation" # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    ExtVolt0 = "volt 0" # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    ExtVolt1 = "volt 1" # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    ExtVolt2 = "volt 2" # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    ExtVolt3 = "volt 3" # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    ExtVolt4 = "volt 4" # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    ExtVolt5 = "volt 5" # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    SBE38temperature = "SBE38 temperature" # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    wetlabs0 = "wetlabs - channel 0" # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    wetlabs1 = "wetlabs - channel 1" # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    wetlabs2 = "wetlabs - channel 2" # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    GTDpressure = "GTD pressure" # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    GTDtemperature = "GTD temperature" # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    GTDpressure2 = "GTD pressure - sensor 2" # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    GTDtemperature2 = "GTD temperature - sensor 2" # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    optodeOxygen = "optode oxygen" # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    SBE63phase = "SBE63 oxygen phase" # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    SBE63temperature = "SBE63 oxygen temperature" # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    dateTime = "date time" # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75


# export type HexDataTypeStrings = keyof typeof HexDataTypes;

HEX_LENGTH = {
    "temperature": 6,
    "conductivity": 6,
    "pressure": 6,
    "temperatureCompensation": 4,
    "voltage": 4,
    "SBE38temperature": 6,
    "wetlabsSingleSensor": 4,  # There are three of these for each WL Sensor
    "GTDpressure": 8,
    "GTDtemperature": 6,
    "optodeOxygen": 6,
    "SBE63phase": 6,
    "SBE63temperature": 6,
    "SeaFETVint": 6,
    "SeaFETVext": 6,
    "SeaOWLChannel": 4,  # There are three of these for each SeaOWL
    "time": 8,
}


class Sensors(Enum):
    """Available sensors to read hex data from"""

    Temperature = "Temperature" # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    Conductivity = "Conductivity" # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    Pressure = "Pressure" # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    ExtVolt0 = "ExtVolt0" # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    ExtVolt1 = "ExtVolt1" # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    ExtVolt2 = "ExtVolt2" # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    ExtVolt3 = "ExtVolt3" # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    ExtVolt4 = "ExtVolt4" # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    ExtVolt5 = "ExtVolt5" # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    WETLABS = "WETLABS" # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    GTD = "GTD" # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    DualGTD = "DualGTD" # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    OPTODE = "OPTODE" # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    SBE63 = "SBE63" # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    SBE38 = "SBE38" # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    SeaFET = "SeaFET" # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75


@dataclass
class MeasurementSeries:
    """Container for measurement data."""

    label: str
    description: str
    units: str
    start_time: date
    values: np.ndarray


@dataclass
class InstrumentData:
    """Container for instrument data parsed from a CNV file."""

    measurements: Dict[str, MeasurementSeries]
    interval_s: float
    latitude: float
    start_time: date
    sample_count: int


def cnv_to_instrument_data(filepath: Path) -> InstrumentData:
    """Import the data from a .cnv file and put it into an
    InstrumentData object.

    :param filepath: the path to the .cnv file to be imported

    :return: the imported data from the .cnv file

    """

    data = InstrumentData(
        measurements={},
        interval_s=None,
        latitude=0.0,
        start_time=None,
        sample_count=None,
    )

    n = 0

    logger.info("Unpacking instrument data from file: %s", filepath)

    with open(filepath, "r", encoding="utf-8") as cnv:
        for line in cnv:
            if line.startswith("*") or line.startswith("#"):
                if line.startswith("# nvalues = "):
                    data.sample_count = int(line[line.find("= ") + 2 : line.find("\n")])
                elif line.startswith("# name "):
                    label = line[line.find("= ") + 2 : line.find(":")]
                    left_bracket = line.find(" [")
                    if left_bracket > 0:
                        description = line[line.find(": ") + 2 : left_bracket]
                        units = line[line.find("[") + 1 : line.find("]")]
                    else:
                        description = line[line.find(": ") + 2 : line.find("\n")].strip()
                        units = ""

                    data.measurements[label] = MeasurementSeries(
                        label=label,
                        description=description,
                        units=units,
                        start_time=None,
                        values=np.zeros(data.sample_count),
                    )

                elif line.startswith("# interval = "):
                    interval = float(
                        line[line.find(": ") + 2 : line.find("\n")]
                    )  # TODO: fix for minutes, hours, etc
                    data.interval_s = interval

                elif line.startswith("# start_time = "):
                    end = min(line.find("\n"), line.find(" ["))
                    date_string = line[line.find("= ") + 2 : end]
                    start_time = datetime.strptime(date_string, "%b %d %Y %H:%M:%S")
                    data.start_time = start_time
                    for measurement in data.measurements.values():
                        measurement.start_time = start_time

                elif line.startswith("** Latitude: "):
                    latitude_parts = line[line.find(": ") + 2 :].split()
                    data.latitude = float(latitude_parts[0]) + float(latitude_parts[1]) / 60.0
                    # TODO: add higher priority latitude to individual measurement series where necessary # pylint: disable=line-too-long

            else:
                values = fix_exponents(" -".join(line.split("-")).split())

                if len(values) > 0:
                    for value, measurement in zip(values, data.measurements.values()):
                        measurement.values[n] = float(value)
                    n += 1
    return data


def fix_exponents(values: List[str]) -> List[str]:
    """Fixes flag values and other numbers with negative exponents.

    The fixes are performed by merging numbers that end in 'e' with the
    following number in the list (the exponent), then deleting the
    exponent

    :param values: List of strings representing numbers

    :return: List of strings where eponents have been fixed
    """

    del_indices = [n + 1 for n, value in enumerate(values) if value.endswith("e")]
    for n in del_indices:
        values[n - 1] = f"{values[n-1]}{values[n]}"
    return np.delete(values, del_indices)


def read_hex_file(
    filepath: str,
    instrument_type: InstrumentType,
    enabled_sensors: List[Sensors],
    moored_mode=False,
) -> pd.DataFrame:
    """Reads a .hex file from a 19plus or 37SM

    :param filepath: path to the .hex file
    :param instrument_type: the instrument that generated the .hex file
    :param enabled_sensors: list of sensors that were enabled on the
        instrument
    :param moored_mode: whether the 19 plus was in moored or profiling
        mode, defaults to False
    :return: a pandas DataFrame with the hex data
    """
    data_count = 0
    is_data = False

    # iterating over file twice in order to preallocate arrays
    file = open(filepath, "r", encoding="utf-8")
    for line in file:
        if is_data and not (line == "" or line.startswith("\n") or line.startswith("\r")):
            data_count += 1
        if line.startswith("*END*"):
            is_data = True

    data_length = data_count
    file.seek(0)
    data = pd.DataFrame()
    data_count = 0
    is_data = False

    for line in file:
        if is_data and not (line == "" or line.startswith("\n") or line.startswith("\r")):
            if data_count == 0:
                data = preallocate_dataframe(
                    instrument_type, line, enabled_sensors, moored_mode, data_length
                )
            hex_data = read_hex(instrument_type, line, enabled_sensors, moored_mode)
            data.iloc[data_count] = pd.DataFrame(hex_data, index=[data_count])
            data_count += 1
        if line.startswith("*END*"):
            is_data = True

    file.close()

    return data


def preallocate_dataframe(
    instrument_type: InstrumentType,
    line: str,
    enabled_sensors: List[Sensors],
    moored_mode: bool,
    data_length: int,
) -> pd.DataFrame:
    """Prefills a pandas DataFrame with zeros for the instrument data

    :param instrument_type: the instrument type
    :param line: TODO: remove in TKIT-63
    :param enabled_sensors: list of sensors that were enabled on the
        instrument
    :param moored_mode: whether the 19 plus was in moored or profiling
        mode
    :param data_length: the number of rows of data in the hex file

    :return: a dataframe fill of zeros
    """
    sensors = {}
    hex_data = read_hex(instrument_type, line, enabled_sensors, moored_mode)
    hex_keys = pd.DataFrame(hex_data, index=[0]).columns
    for key in hex_keys:
        sensors[key] = np.zeros(data_length)

    return pd.DataFrame(sensors)


def read_hex(
    instrument_type: InstrumentType,
    hex: str,
    enabled_sensors: List[Sensors],
    moored_mode=False,
) -> dict:
    """Converts an instrument data hex string into engineering units.

    :param instrument_type: determines how units are converted
    :param hex: one line from a hex data file
    :param enabled_sensors: mooredMode parses time for 19plus in moored
        mode if true
    :param moored_mode: array of Sensors that are enabled. For 37 this
        is always temperature, conductivity, pressure. Defaults to False

    :return: the sensor values in engineering units that were extracted
        from the input hex string
    """

    if instrument_type == InstrumentType.SBE19Plus:
        return read_SBE19plus_format_0(hex, enabled_sensors, moored_mode)

    if instrument_type == InstrumentType.SBE37SM:
        return read_SBE37SM_format_0(hex, enabled_sensors)


def read_SBE19plus_format_0(hex: str, enabled_sensors: List[Sensors], moored_mode=False) -> dict: # pylint: disable=invalid-name, line-too-long #TODO: change this to be snake_case for TKIT-75
    """Converts a 19plus V2 data hex string into engineering units.

    :param hex: one line from a hex data file
    :param enabled_sensors: array of Sensors that are enabled. For 37
        this is always temperature, conductivity, pressure. Defaults to
        False
    :param moored_mode: parses time for 19plus in moored mode if true

    :return: the 19plus V2 sensor values in engineering units that were
            extracted from the input hex string

    :raises RuntimeWarning: if the hex string length does not match the
        expected length
    """

    results = {}
    n = 0
    for sensor in Sensors:
        if sensor in enabled_sensors:
            if sensor == Sensors.Temperature:
                results[HexDataTypes.temperature.value] = int(
                    hex[n : HEX_LENGTH["temperature"]], 16
                )
                n += HEX_LENGTH["temperature"]

            if sensor == Sensors.Conductivity:
                results[HexDataTypes.conductivity.value] = (
                    int(hex[n : n + HEX_LENGTH["conductivity"]], 16) / 256
                )
                n += HEX_LENGTH["conductivity"]

            if sensor == Sensors.Pressure:  # TODO: add conversion for quartz pressure sensors
                results[HexDataTypes.pressure.value] = int(hex[n : n + HEX_LENGTH["pressure"]], 16)
                n += HEX_LENGTH["pressure"]
                result = (
                    int(hex[n : n + HEX_LENGTH["temperatureCompensation"]], 16) / COUNTS_TO_VOLTS
                )
                results[HexDataTypes.temperatureCompensation.value] = result
                n += HEX_LENGTH["temperatureCompensation"]

            if sensor in [
                Sensors.ExtVolt0,
                Sensors.ExtVolt1,
                Sensors.ExtVolt2,
                Sensors.ExtVolt3,
                Sensors.ExtVolt4,
                Sensors.ExtVolt5,
            ]:
                result = int(hex[n : n + HEX_LENGTH["voltage"]], 16) / COUNTS_TO_VOLTS
                results[HexDataTypes[Sensors[sensor.value].value].value] = result
                n += HEX_LENGTH["voltage"]

            if sensor == Sensors.SBE38:
                results[HexDataTypes.SBE38temperature.value] = (
                    int(hex[n : n + HEX_LENGTH["SBE38temperature"]], 16) / 100000 - 10
                )
                n += HEX_LENGTH["SBE38temperature"]

            if sensor == Sensors.WETLABS:
                results[HexDataTypes.wetlabs0.value] = int(
                    hex[n : n + HEX_LENGTH["wetlabsSingleSensor"]], 16
                )
                n += HEX_LENGTH["wetlabsSingleSensor"]

                results[HexDataTypes.wetlabs1.value] = int(
                    hex[n : n + HEX_LENGTH["wetlabsSingleSensor"]], 16
                )
                n += HEX_LENGTH["wetlabsSingleSensor"]

                results[HexDataTypes.wetlabs2.value] = int(
                    hex[n : n + HEX_LENGTH["wetlabsSingleSensor"]], 16
                )
                n += HEX_LENGTH["wetlabsSingleSensor"]

            if sensor == Sensors.GTD:
                results[HexDataTypes.GTDpressure.value] = (
                    int(hex[n : n + HEX_LENGTH["GTDpressure"]], 16) / 10000
                )
                n += HEX_LENGTH["GTDpressure"]
                results[HexDataTypes.GTDtemperature.value] = (
                    int(hex[n : n + HEX_LENGTH["GTDtemperature"]], 16) / 10000 - 10
                )
                n += HEX_LENGTH["GTDtemperature"]

            if sensor == Sensors.DualGTD:
                results[HexDataTypes.GTDpressure.value] = (
                    int(hex[n : n + HEX_LENGTH["GTDpressure"]], 16) / 10000
                )
                n += HEX_LENGTH["GTDpressure"]
                results[HexDataTypes.GTDtemperature.value] = (
                    int(hex[n : n + HEX_LENGTH["GTDtemperature"]], 16) / 10000 - 10
                )
                n += HEX_LENGTH["GTDtemperature"]
                results[HexDataTypes.GTDpressure2.value] = (
                    int(hex[n : n + HEX_LENGTH["GTDpressure"]], 16) / 10000
                )
                n += HEX_LENGTH["GTDpressure"]
                results[HexDataTypes.GTDtemperature2.value] = (
                    int(hex[n : n + HEX_LENGTH["GTDtemperature"]], 16) / 10000 - 10
                )
                n += HEX_LENGTH["GTDtemperature"]

            if sensor == Sensors.OPTODE:
                results[HexDataTypes.optodeOxygen.value] = (
                    int(hex[n : n + HEX_LENGTH["optodeOxygen"]], 16) / 10000 - 10
                )
                n += HEX_LENGTH["optodeOxygen"]

            if sensor == Sensors.SBE63:
                results[HexDataTypes.SBE63phase.value] = (
                    int(hex[n : n + HEX_LENGTH["SBE63phase"]], 16) / 100000 - 10
                )
                n += HEX_LENGTH["SBE63phase"]
                results[HexDataTypes.SBE63temperature.value] = (
                    int(hex[n : n + HEX_LENGTH["SBE63temperature"]], 16) / 1000000 - 1
                )
                n += HEX_LENGTH["SBE63temperature"]

    if moored_mode:
        seconds_since_2000 = int(hex[n : n + HEX_LENGTH["time"]], 16)
        results[HexDataTypes.dateTime.value] = datetime.fromtimestamp(
            seconds_since_2000 + SECONDS_BETWEEN_EPOCH_AND_2000
        )
        n += HEX_LENGTH["time"]

    # Validate hex length. Ensure length matches what is expected based
    # on enabled sensors and moored mode.
    if n != len(hex.split("\n")[0]):
        raise RuntimeWarning(
            "Hex string length does not match expectation based on enabled sensors and moored mode"
        )

    return results


def read_SBE37SM_format_0(hex: str, enabled_sensors: List[Sensors]) -> dict: # pylint: disable=invalid-name, line-too-long #TODO: change this to be snake_case for TKIT-75
    """Converts a 37 family data hex string into engineering units.

    :param hex: one line from a hex data file
    :param enabled_sensors: array of Sensors that are enabled. For 37
        this is always temperature, conductivity, pressure. Defaults to
        False

    :return: the 37 family sensor values in engineering units that were
        extracted from the input hex string
    """
    results = {}
    n = 0
    results[HexDataTypes.temperature.value] = int(hex[n : HEX_LENGTH["temperature"]], 16)
    n += HEX_LENGTH["temperature"]

    results[HexDataTypes.conductivity.value] = (
        int(hex[n : n + HEX_LENGTH["conductivity"]], 16) / 256
    )
    n += HEX_LENGTH["conductivity"]

    if Sensors.SBE63 in enabled_sensors:
        results[HexDataTypes.SBE63phase.value] = (
            int(hex[n : n + HEX_LENGTH["SBE63phase"]], 16) / 100000 - 10
        )
        n += HEX_LENGTH["SBE63phase"]
        results[HexDataTypes.SBE63temperature.value] = (
            int(hex[n : n + HEX_LENGTH["SBE63temperature"]], 16) / 1000000 - 1
        )
        n += HEX_LENGTH["SBE63temperature"]

    if Sensors.Pressure in enabled_sensors:
        results[HexDataTypes.pressure.value] = int(hex[n : n + HEX_LENGTH["pressure"]], 16)
        n += HEX_LENGTH["pressure"]
        result = int(hex[n : n + HEX_LENGTH["temperatureCompensation"]], 16)
        results[HexDataTypes.temperatureCompensation.value] = result
        n += HEX_LENGTH["temperatureCompensation"]

    seconds_since_2000 = int(hex[n : n + HEX_LENGTH["time"]], 16)
    results[HexDataTypes.dateTime.value] = datetime.fromtimestamp(
        seconds_since_2000 + SECONDS_BETWEEN_EPOCH_AND_2000
    )
    n += HEX_LENGTH["time"]

    return results
