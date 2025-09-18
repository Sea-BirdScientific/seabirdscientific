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
#   read_SBE19plus_format_0
#   read_SBE37SM_format_0
#   read_SBE39plus_format_0
#   read_seafet_format_0
#   read_SBE911plus_format_0

# Native imports
import builtins
from enum import Enum
from dataclasses import dataclass
from datetime import date, datetime, timedelta
from logging import getLogger
from typing import List, Dict, Optional, Union
from pathlib import Path
import warnings

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
    SBE37SI = "37-SI"
    SBE37SIP = "37-SIP"
    SBE37IM = "37-IM"
    SBE37IMP = "37-IMP"
    SBE37IMPODO = "37-IMP-ODO"
    SBE19Plus = "19plus"  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    SBE16Plus = "16plus"
    SBE39Plus = "39plus"
    SBE39PlusIM = "39plus-IM"
    SBE911Plus = "911plus"
    SeaFET2 = "SeaFET2"
    SeapHox2 = "SeapHox2"
    HydroCAT = "HydroCAT"
    HydroCATODO = "HydroCAT-ODO"


class HexDataTypes(Enum):
    """Possible data types in hex files"""

    temperature = (  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
        "temperature"
    )
    secondaryTemperature = (  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
        "secondary temperature"
    )
    conductivity = (  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
        "conductivity"
    )
    secondaryConductivity = (  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
        "secondary conductivity"
    )

    pressure = "pressure"  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    digiquartzPressure = "digiquartz pressure"  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    temperatureCompensation = "temperature compensation"  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    ExtVolt0 = "volt 0"  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    ExtVolt1 = "volt 1"  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    ExtVolt2 = "volt 2"  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    ExtVolt3 = "volt 3"  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    ExtVolt4 = "volt 4"  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    ExtVolt5 = "volt 5"  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    ExtVolt6 = "volt 6"
    ExtVolt7 = "volt 7"
    SBE38temperature = "SBE38 temperature"  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    wetlabs0 = "wetlabs - channel 0"  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    wetlabs1 = "wetlabs - channel 1"  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    wetlabs2 = "wetlabs - channel 2"  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    GTDpressure = (  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
        "GTD pressure"
    )
    GTDtemperature = (  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
        "GTD temperature"
    )
    GTDpressure2 = "GTD pressure - sensor 2"  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    GTDtemperature2 = "GTD temperature - sensor 2"  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    optodeOxygen = (  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
        "optode oxygen"
    )
    SBE63phase = "SBE63 oxygen phase"  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    SBE63temperature = "SBE63 oxygen temperature"  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    dateTime = "date time"  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    # NMEA Devices
    nmeaTime = (  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
        "NMEA Date Time"
    )
    nmeaLatitude = (  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
        "NMEA Latitude"
    )
    nmeaLongitude = (  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
        "NMEA Longitude"
    )
    nmeaDepth = "nmea depth"
    statusAndSign = (  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
        "status and sign"
    )
    vrsInternal = (  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
        "vrs internal"
    )
    vrsExternal = (  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
        "vrs external"
    )
    pHtemperature = (  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
        "ph temperature"
    )
    vk = "vk"  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    ib = "ib"  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    ik = "ik"  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    relativeHumidity = (  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
        "relative humidity"
    )
    internalTemperature = (  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
        "internal temperature"
    )
    errorFlag = (  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
        "error flag"
    )
    surfacePAR = "surface par"
    SBE911PumpStatus = "SBE911 pump status"
    SBE911BottomContactStatus = "SBE911 bottom contact status"
    SBE911ConfirmStatus = "SBE911 confirm status"
    SBE911ModemStatus = "SBE911 modem status"
    dataIntegrity = "data integrity"
    systemTime = "system time"


# export type HexDataTypeStrings = keyof typeof HexDataTypes;

HEX_LENGTH = {
    "temperature": 6,
    "conductivity": 6,
    "pressure": 6,
    "temperatureCompensation": 4,
    "voltage": 4,
    "SBE911Voltage": 3,
    "SBE911SPAR": 3,
    "SBE911TemperatureCompensation": 3,
    "SBE911Status": 1,
    "SBE911DataIntegrity": 2,
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
    "vrsExternal": 6,
    "vrsInternal": 6,
    "pHtemperature": 6,
    "vk": 6,
    "ib": 6,
    "ik": 6,
    "relativeHumidity": 3,
    "internalTemperature": 3,
    "errorFlag": 4,
    # nmea devices
    "nmeaLatitude": 6,
    "nmeaLongitude": 6,
    "nmeaTime": 8,
    "nmeaLocation": 14,
    "statusAndSign": 2,
    "systemTime": 8,
}


class Sensors(Enum):
    """Available sensors to read hex data from"""

    Temperature = (  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
        "Temperature"
    )
    SecondaryTemperature = "SecondaryTemperature"  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    Conductivity = (  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
        "Conductivity"
    )
    SecondaryConductivity = "SecondaryConductivity"  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    Pressure = "Pressure"  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    ExtVolt0 = "ExtVolt0"  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    ExtVolt1 = "ExtVolt1"  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    ExtVolt2 = "ExtVolt2"  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    ExtVolt3 = "ExtVolt3"  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    ExtVolt4 = "ExtVolt4"  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    ExtVolt5 = "ExtVolt5"  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    ExtVolt6 = "ExtVolt6"  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    ExtVolt7 = "ExtVolt7"
    WETLABS = "WETLABS"  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    GTD = "GTD"  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    DualGTD = "DualGTD"  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    OPTODE = "OPTODE"  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    SBE63 = "SBE63"  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    SBE38 = "SBE38"  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    SeaFET = "SeaFET"  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    SPAR = "SPAR"  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    # nmea devices
    nmeaLatitude = (  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
        "nmeaLatitude"
    )
    nmeaLongitude = (  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
        "nmeaLongitude"
    )
    statusAndSign = (  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
        "StatusAndSign"
    )
    nmeaTime = "nmeaTime"  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    nmeaLocation = (
        "nmeaLocation"  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    )
    nmeaDepth = (
        "nmeaDepth"  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    )
    SystemTime = (
        "systemTime"  # pylint: disable=invalid-name # change enums to UPPER_CASE for TKIT-75
    )


@dataclass
class MeasurementSeries:
    """Container for measurement data."""

    label: str
    description: str
    units: str
    start_time: Optional[date]
    values: np.ndarray


@dataclass
class InstrumentData:
    """Container for instrument data parsed from a CNV file."""

    measurements: Dict[str, MeasurementSeries]
    """
    Dictionary of MeasurementSeries by their label.

    Note: duplicate labels can occur and cnv_to_instrument_data will numerically
    increment the duplicate labels, e.g. the second "depSM" becomes "depSM1".
    """
    interval_s: Optional[float]
    latitude: float
    start_time: Optional[date]
    sample_count: Optional[int]

    def _to_dataframe(self):
        measurements = {k: v.values for k, v in self.measurements.items()}
        return pd.DataFrame(measurements)


# pylint: disable=too-many-branches # TODO: Fix this
def cnv_to_instrument_data(filepath: Path) -> InstrumentData:
    """
    Import the data from a .cnv file and put it into an InstrumentData object.

    Duplicate labels will be incremented for InstrumentData.measurements keys.
    For example, the second "depSM" becomes "depSM1".
    However, the MeasurementSeries.label will be the original label.

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

                    num_values = data.sample_count or 0  # num_values to 0 if sample_count is None

                    key = label
                    # check for label collision
                    # if collision, increment until find available key
                    i = 1
                    while key in data.measurements:
                        key = f"{label}_{i}"
                        i += 1

                    data.measurements[key] = MeasurementSeries(
                        label=label,
                        description=description,
                        units=units,
                        start_time=None,
                        values=np.zeros(num_values),
                    )

                    if key != label:
                        logger.warning(
                            'duplicate measurement "%s" will use key "%s" in measurements dict',
                            label,
                            key,
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
                    # TODO: add higher priority latitude to individual measurement series
                    # where necessary

            else:
                values = fix_exponents(" -".join(line.split("-")).split())

                if len(values) > 0:
                    for value, measurement in zip(values, data.measurements.values()):
                        measurement.values[n] = float(value)
                    n += 1
    return data


def fix_exponents(values: List[str]) -> List[str]:
    """Fixes flag values and other numbers with negative exponents.
    This is necessary because sometimes there is only a minus sign
    separating two values in a cnv file. So we split values on the minus
    sign which also splits negative exponents (e.g. 1e-2 becomes 1e, -2).
    This function repairs the exponents by merging numbers that end in
    'e' with the following number in the list (e.g. 1e, -2 becomes 1e-2),
    then removes the extra exponent from the list.

    :param values: List of strings representing numbers

    :return: List of strings where eponents have been fixed
    """

    del_indices = [n + 1 for n, value in enumerate(values) if value.endswith("e")]
    for n in del_indices:
        values[n - 1] = f"{values[n-1]}{values[n]}"
    new_values = list(np.delete(values, del_indices))
    return new_values


def read_hex_file(
    filepath: str,
    instrument_type: InstrumentType,
    enabled_sensors: List[Sensors],
    moored_mode=False,
    is_shallow=True,
) -> pd.DataFrame:
    """Reads a .hex file

    :param filepath: path to the .hex file
    :param instrument_type: the instrument that generated the .hex file
    :param enabled_sensors: list of sensors that were enabled on the
        instrument
    :param moored_mode: whether the instrument was in moored or profiling
        mode, defaults to False
    :return: a pandas DataFrame with the hex data
    """
    data_count = 0
    is_data = False

    # iterating over file twice in order to preallocate arrays
    # pylint: disable=consider-using-with # TODO: Fix this
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
            hex_data = read_hex(instrument_type, line, enabled_sensors, moored_mode, is_shallow)
            for key, value in hex_data.items():
                data.loc[data_count, key] = value
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
    # hex_keys = pd.DataFrame(hex_data, index=[0]).columns
    for key, value in hex_data.items():
        if isinstance(value, datetime):
            sensors[key] = pd.date_range(start="2000-01-01", end="2000-01-02", periods=data_length)
        else:
            sensors[key] = np.zeros(data_length)

    return pd.DataFrame(sensors)


def read_hex(
    instrument_type: InstrumentType,
    hex_segment: str = "",
    enabled_sensors: List[Sensors] | None = None,
    moored_mode=False,
    is_shallow=True,
    frequency_channels_supressed=0,
    voltage_words_supressed=0,
    hex=hex,  # pylint: disable=redefined-builtin
) -> dict:
    """Converts an instrument data hex string into engineering units.

    :param instrument_type: determines how units are converted
    :param hex_segment: one line from a hex data file
    :param enabled_sensors: mooredMode parses time for 19plus in moored
        mode if true
    :param moored_mode: array of Sensors that are enabled. For 37 this
        is always temperature, conductivity, pressure. Defaults to False
    :param hex: Deprecated, use hex_segment

    :return: the sensor values in engineering units that were extracted
        from the input hex string
    """
    if hex is not builtins.hex:
        warnings.warn("hex is deprecated, use hex_segment", DeprecationWarning)

    match instrument_type:
        case InstrumentType.SBE19Plus:
            return read_SBE19plus_format_0(hex_segment, enabled_sensors, moored_mode)
        case InstrumentType.SBE16Plus:
            return read_SBE19plus_format_0(hex_segment, enabled_sensors, True)

        case (
            InstrumentType.SBE37IM
            | InstrumentType.SBE37IMP
            | InstrumentType.SBE37IMPODO
            | InstrumentType.SBE37IM
            | InstrumentType.SBE37SI
            | InstrumentType.SBE37SIP
            | InstrumentType.SBE37SM
            | InstrumentType.SBE37SMP
            | InstrumentType.SBE37SMPODO
            | InstrumentType.HydroCAT
            | InstrumentType.HydroCATODO
        ):
            return read_SBE37SM_format_0(hex_segment, enabled_sensors)

        case InstrumentType.SBE39Plus | InstrumentType.SBE39PlusIM:
            return read_SBE39plus_format_0(hex_segment, enabled_sensors)

        case InstrumentType.SeaFET2 | InstrumentType.SeapHox2:
            return read_seafet_format_0(hex_segment, instrument_type, is_shallow)
        case InstrumentType.SBE911Plus:
            return read_SBE911plus_format_0(
                hex_segment, enabled_sensors, frequency_channels_supressed, voltage_words_supressed
            )

        case _:
            raise ValueError(f"Unsupported instrument type: {instrument_type}")

    return {}


def read_SBE39plus_format_0(
    hex_segment: str = "",
    enabled_sensors: List[Sensors] | None = None,
    hex=builtins.hex,
) -> Dict[str, Union[int, float, datetime]]:
    """Converts a 39plus data hex string into engineering units.

    :param hex_segment: one line from a hex data file
    :param enabled_sensors: array of Sensors that are enabled
    :param hex: Deprecated, use hex_segment

    :return: the 39plus sensor values in engineering units that were
        extracted from the input hex string
    """
    if hex is not builtins.hex:
        warnings.warn("hex is deprecated, use hex_segment", DeprecationWarning)

    results: Dict[str, Union[int, float, datetime]] = {}
    n = 0

    # Datetime (naive, no timezone conversion)
    seconds_since_2000 = int(hex_segment[n : n + HEX_LENGTH["time"]], 16)
    results[HexDataTypes.dateTime.value] = datetime(1970, 1, 1) + timedelta(
        seconds=seconds_since_2000 + SECONDS_BETWEEN_EPOCH_AND_2000
    )
    n += HEX_LENGTH["time"]

    # Temperature
    results[HexDataTypes.temperature.value] = int(
        hex_segment[n : n + HEX_LENGTH["temperature"]], 16
    )
    n += HEX_LENGTH["temperature"]

    # Pressure and temperature compensation
    if enabled_sensors and Sensors.Pressure in enabled_sensors:
        results[HexDataTypes.pressure.value] = int(hex_segment[n : n + HEX_LENGTH["pressure"]], 16)
        n += HEX_LENGTH["pressure"]

        temp_comp = int(hex_segment[n : n + HEX_LENGTH["temperatureCompensation"]], 16)
        results[HexDataTypes.temperatureCompensation.value] = temp_comp
        n += HEX_LENGTH["temperatureCompensation"]

    # Validate hex length
    if n != len(hex_segment.strip()):
        raise ValueError("Hex string length does not match expectation based on enabled sensors")

    return results


def read_seafet_format_0(
    hex_segment: str = "",
    instrument_type: InstrumentType = None,
    is_shallow: bool = True,
    hex=builtins.hex,
) -> Dict[str, Union[int, float, datetime]]:
    """Converts a SeaFET2 or SeapHox2 hex string into engineering units.

    :param hex_segment: one line from a hex data file
    :param instrument_type: InstrumentType.SeaFET2 or InstrumentType.SeapHox2
    :param is_shallow: if True, include internal pH and pH reference temperature
    :param hex: Deprecated, use hex_segment

    :return: sensor values in engineering units extracted from the hex string
    """
    if hex is not builtins.hex:
        warnings.warn("hex is deprecated, use hex_segment", DeprecationWarning)

    if instrument_type not in (InstrumentType.SeaFET2, InstrumentType.SeapHox2):
        raise ValueError(
            f"In read_seafet_format_0 {instrument_type} is not recognized as a SeaFET2 or SeapHox2"
        )

    results: Dict[str, Union[int, float, datetime]] = {}
    n = 0

    # SeapHox2 specific values
    if instrument_type == InstrumentType.SeapHox2:
        results[HexDataTypes.temperature.value] = int(
            hex_segment[n : n + HEX_LENGTH["temperature"]], 16
        )
        n += HEX_LENGTH["temperature"]

        results[HexDataTypes.conductivity.value] = (
            int(hex_segment[n : n + HEX_LENGTH["conductivity"]], 16) / 256
        )
        n += HEX_LENGTH["conductivity"]

        results[HexDataTypes.pressure.value] = int(hex_segment[n : n + HEX_LENGTH["pressure"]], 16)
        n += HEX_LENGTH["pressure"]

        results[HexDataTypes.temperatureCompensation.value] = int(
            hex_segment[n : n + HEX_LENGTH["temperatureCompensation"]], 16
        )
        n += HEX_LENGTH["temperatureCompensation"]

        results[HexDataTypes.SBE63phase.value] = (
            int(hex_segment[n : n + HEX_LENGTH["SBE63phase"]], 16) / 100000 - 10
        )
        n += HEX_LENGTH["SBE63phase"]

        results[HexDataTypes.SBE63temperature.value] = (
            int(hex_segment[n : n + HEX_LENGTH["SBE63temperature"]], 16) / 1000000 - 1
        )
        n += HEX_LENGTH["SBE63temperature"]

    # External pH
    results[HexDataTypes.vrsExternal.value] = int(
        hex_segment[n : n + HEX_LENGTH["vrsExternal"]], 16
    )
    n += HEX_LENGTH["vrsExternal"]

    if is_shallow:
        # Internal pH
        results[HexDataTypes.vrsInternal.value] = int(
            hex_segment[n : n + HEX_LENGTH["vrsInternal"]], 16
        )
        n += HEX_LENGTH["vrsInternal"]

        # pH reference temperature
        results[HexDataTypes.pHtemperature.value] = int(
            hex_segment[n : n + HEX_LENGTH["pHtemperature"]], 16
        )
        n += HEX_LENGTH["pHtemperature"]

    # Other sensors
    results[HexDataTypes.vk.value] = int(hex_segment[n : n + HEX_LENGTH["vk"]], 16)
    n += HEX_LENGTH["vk"]

    results[HexDataTypes.ib.value] = int(hex_segment[n : n + HEX_LENGTH["ib"]], 16)
    n += HEX_LENGTH["ib"]

    results[HexDataTypes.ik.value] = int(hex_segment[n : n + HEX_LENGTH["ik"]], 16)
    n += HEX_LENGTH["ik"]

    results[HexDataTypes.relativeHumidity.value] = int(
        hex_segment[n : n + HEX_LENGTH["relativeHumidity"]] + "0", 16
    )
    n += HEX_LENGTH["relativeHumidity"]

    results[HexDataTypes.internalTemperature.value] = int(
        hex_segment[n : n + HEX_LENGTH["internalTemperature"]] + "0", 16
    )
    n += HEX_LENGTH["internalTemperature"]

    # Datetime (naive, no timezone conversion)
    seconds_since_2000 = int(hex_segment[n : n + HEX_LENGTH["time"]], 16)
    results[HexDataTypes.dateTime.value] = datetime(2000, 1, 1) + timedelta(
        seconds=seconds_since_2000
    )
    n += HEX_LENGTH["time"]

    # Error flag
    results[HexDataTypes.errorFlag.value] = int(hex_segment[n : n + HEX_LENGTH["errorFlag"]], 16)
    n += HEX_LENGTH["errorFlag"]

    # Validate hex length
    if n != len(hex_segment.strip()):
        raise ValueError(
            "Hex string length does not match expectation based on instrument_type and is_shallow"
        )

    return results


# TODO: change this to be snake_case for TKIT-75
def read_SBE911plus_format_0(  # pylint: disable=invalid-name
    hex_segment: str = "",
    enabled_sensors: list["Sensors"] | None = None,
    frequency_channels_suppressed: int = 0,
    voltage_words_suppressed: int = 0,
) -> dict[str, Union[int, float, datetime]]:
    """Converts a 911Plus hex string into engineering units.

    :param hex_segment: one line from a hex data file
    :param enabled_sensors: list of enabled Sensors
    :param frequency_channels_suppressed: number of suppressed frequency channels
    :param voltage_words_suppressed: number of suppressed voltage words

    :return: dictionary of sensor values in engineering units
    """
    if enabled_sensors is None:
        enabled_sensors = []

    results: dict[str, Union[int, float, datetime]] = {}
    n = 0

    # Temperature
    results[HexDataTypes.temperature.value] = frequency_from_3_bytes(
        hex_segment[n : n + HEX_LENGTH["temperature"]]
    )
    n += HEX_LENGTH["temperature"]

    # Conductivity
    results[HexDataTypes.conductivity.value] = frequency_from_3_bytes(
        hex_segment[n : n + HEX_LENGTH["conductivity"]]
    )
    n += HEX_LENGTH["conductivity"]

    # Digiquartz Pressure
    results[HexDataTypes.digiquartzPressure.value] = frequency_from_3_bytes(
        hex_segment[n : n + HEX_LENGTH["pressure"]]
    )
    n += HEX_LENGTH["pressure"]

    # Secondary temperature
    if Sensors.SecondaryTemperature in enabled_sensors:
        results[HexDataTypes.secondaryTemperature.value] = frequency_from_3_bytes(
            hex_segment[n : n + HEX_LENGTH["temperature"]]
        )
        n += HEX_LENGTH["temperature"]
    elif frequency_channels_suppressed <= 1:
        n += HEX_LENGTH["temperature"]

    # Secondary conductivity
    if Sensors.SecondaryConductivity in enabled_sensors:
        results[HexDataTypes.secondaryConductivity.value] = frequency_from_3_bytes(
            hex_segment[n : n + HEX_LENGTH["conductivity"]]
        )
        n += HEX_LENGTH["conductivity"]
    elif frequency_channels_suppressed == 0:
        n += HEX_LENGTH["conductivity"]

    # Voltage channels (suppressed in pairs)
    if Sensors.ExtVolt0 in enabled_sensors or Sensors.ExtVolt1 in enabled_sensors:
        results[HexDataTypes.ExtVolt0.value], results[HexDataTypes.ExtVolt1.value] = (
            voltages_from_3_bytes(hex_segment[n : n + HEX_LENGTH["SBE911Voltage"] * 2])
        )
        n += HEX_LENGTH["SBE911Voltage"] * 2
    elif voltage_words_suppressed <= 3:
        # NotInUse volt channel that was not suppressed
        n += HEX_LENGTH["SBE911Voltage"] * 2

    if Sensors.ExtVolt2 in enabled_sensors or Sensors.ExtVolt3 in enabled_sensors:
        results[HexDataTypes.ExtVolt2.value], results[HexDataTypes.ExtVolt3.value] = (
            voltages_from_3_bytes(hex_segment[n : n + HEX_LENGTH["SBE911Voltage"] * 2])
        )
        n += HEX_LENGTH["SBE911Voltage"] * 2
    elif voltage_words_suppressed <= 2:
        # NotInUse volt channel that was not suppressed
        n += HEX_LENGTH["SBE911Voltage"] * 2

    if Sensors.ExtVolt4 in enabled_sensors or Sensors.ExtVolt5 in enabled_sensors:
        results[HexDataTypes.ExtVolt4.value], results[HexDataTypes.ExtVolt5.value] = (
            voltages_from_3_bytes(hex_segment[n : n + HEX_LENGTH["SBE911Voltage"] * 2])
        )
        n += HEX_LENGTH["SBE911Voltage"] * 2
    elif voltage_words_suppressed <= 1:
        # NotInUse volt channel that was not suppressed
        n += HEX_LENGTH["SBE911Voltage"] * 2

    if Sensors.ExtVolt6 in enabled_sensors or Sensors.ExtVolt7 in enabled_sensors:
        results[HexDataTypes.ExtVolt6.value], results[HexDataTypes.ExtVolt7.value] = (
            voltages_from_3_bytes(hex_segment[n : n + HEX_LENGTH["SBE911Voltage"] * 2])
        )
        n += HEX_LENGTH["SBE911Voltage"] * 2
    elif voltage_words_suppressed == 0:
        # NotInUse volt channel that was not suppressed
        n += HEX_LENGTH["SBE911Voltage"] * 2

    # Surface PAR
    if Sensors.SPAR in enabled_sensors:
        n += 3  # unused bits
        results[HexDataTypes.surfacePAR.value] = (
            int(hex_segment[n : n + HEX_LENGTH["SBE911SPAR"]], 16) / 819
        )
        n += HEX_LENGTH["SBE911SPAR"]

    # NMEA Location
    if Sensors.nmeaLocation in enabled_sensors:
        hex_loc = hex_segment[n : n + HEX_LENGTH["nmeaLocation"]]

        def get_lat_lon(hex_segment: str) -> tuple[float, float]:
            byte0 = int(hex_segment[0:2], 16)
            byte1 = int(hex_segment[2:4], 16)
            byte2 = int(hex_segment[4:6], 16)
            byte3 = int(hex_segment[6:8], 16)
            byte4 = int(hex_segment[8:10], 16)
            byte5 = int(hex_segment[10:12], 16)
            byte6 = format(int(hex_segment[12:14], 16), "08b")
            lat_sign = 1 if byte6[0] == "0" else -1
            lon_sign = 1 if byte6[1] == "0" else -1
            lat = lat_sign * (byte0 * 65536 + byte1 * 256 + byte2) / 50000
            lon = lon_sign * (byte3 * 65536 + byte4 * 256 + byte5) / 50000
            return lat, lon

        results[HexDataTypes.nmeaLatitude.value], results[HexDataTypes.nmeaLongitude.value] = (
            get_lat_lon(hex_loc)
        )
        n += HEX_LENGTH["nmeaLocation"]

    # NMEA Depth
    if Sensors.nmeaDepth in enabled_sensors:
        results[HexDataTypes.nmeaDepth.value] = int(
            hex_segment[n : n + HEX_LENGTH["nmeaDepth"]], 16
        )
        n += HEX_LENGTH["nmeaDepth"]

    # NMEA Time
    if Sensors.nmeaTime in enabled_sensors:
        seconds_since_2000 = int(
            reverse_hex_bytes(hex_segment[n : n + HEX_LENGTH["nmeaTime"]]), 16
        )
        results[HexDataTypes.nmeaTime.value] = datetime(2000, 1, 1) + timedelta(
            seconds=seconds_since_2000
        )
        n += HEX_LENGTH["nmeaTime"]

    # Temperature compensation
    results[HexDataTypes.temperatureCompensation.value] = int(
        hex_segment[n : n + HEX_LENGTH["SBE911TemperatureCompensation"]], 16
    )
    n += HEX_LENGTH["SBE911TemperatureCompensation"]

    # Status bits
    status_bin = format(int(hex_segment[n : n + HEX_LENGTH["SBE911Status"]], 16), "04b")
    results[HexDataTypes.SBE911PumpStatus.value] = int(status_bin[0])
    results[HexDataTypes.SBE911BottomContactStatus.value] = int(status_bin[1])
    results[HexDataTypes.SBE911ConfirmStatus.value] = int(status_bin[2])
    results[HexDataTypes.SBE911ModemStatus.value] = int(status_bin[3])
    n += HEX_LENGTH["SBE911Status"]

    # Data integrity
    results[HexDataTypes.dataIntegrity.value] = int(
        hex_segment[n : n + HEX_LENGTH["SBE911DataIntegrity"]], 16
    )
    n += HEX_LENGTH["SBE911DataIntegrity"]

    # System time
    if Sensors.SystemTime in enabled_sensors:
        seconds_since_1970 = int(
            reverse_hex_bytes(hex_segment[n : n + HEX_LENGTH["systemTime"]]), 16
        )
        results[HexDataTypes.systemTime.value] = datetime(1970, 1, 1) + timedelta(
            seconds=seconds_since_1970
        )
        n += HEX_LENGTH["systemTime"]

    # Validate hex length
    if n != len(hex_segment.strip()):
        raise ValueError("Hex string length does not match expectation based on enabled sensors")

    return results


# TODO: change the following fn name to be snake_case for TKIT-75
# pylint: disable=invalid-name,too-many-branches,too-many-statements # TODO: Fix these
def read_SBE19plus_format_0(
    hex_segment: str = "",
    enabled_sensors: List[Sensors] | None = None,
    moored_mode=False,
    hex=hex,  # pylint: disable=redefined-builtin
) -> Dict[str, Union[float, datetime]]:
    """Converts a 19plus V2 data hex string into engineering units.

    :param hex_segment: one line from a hex data file
    :param enabled_sensors: array of Sensors that are enabled. For 37
        this is always temperature, conductivity, pressure. Defaults to
        False
    :param moored_mode: parses time for 19plus in moored mode if true
    :param hex: Deprecated, use hex_segment

    :return: the 19plus V2 sensor values in engineering units that were
            extracted from the input hex string

    :raises RuntimeWarning: if the hex string length does not match the
        expected length
    """
    if hex is not builtins.hex:
        warnings.warn("hex is deprecated, use hex_segment", DeprecationWarning)

    results: Dict[str, Union[int, float, datetime]] = {}
    n = 0
    for sensor in Sensors:
        if enabled_sensors and sensor in enabled_sensors:
            if sensor == Sensors.Temperature:
                results[HexDataTypes.temperature.value] = int(
                    hex_segment[n : HEX_LENGTH["temperature"]], 16
                )
                n += HEX_LENGTH["temperature"]

            if sensor == Sensors.Conductivity:
                results[HexDataTypes.conductivity.value] = (
                    int(hex_segment[n : n + HEX_LENGTH["conductivity"]], 16) / 256
                )
                n += HEX_LENGTH["conductivity"]

            if sensor == Sensors.Pressure:  # TODO: add conversion for quartz pressure sensors
                results[HexDataTypes.pressure.value] = int(
                    hex_segment[n : n + HEX_LENGTH["pressure"]], 16
                )
                n += HEX_LENGTH["pressure"]
                result = (
                    int(hex_segment[n : n + HEX_LENGTH["temperatureCompensation"]], 16)
                    / COUNTS_TO_VOLTS
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
                result = int(hex_segment[n : n + HEX_LENGTH["voltage"]], 16) / COUNTS_TO_VOLTS
                results[HexDataTypes[Sensors[sensor.value].value].value] = result
                n += HEX_LENGTH["voltage"]

            if sensor == Sensors.SBE38:
                results[HexDataTypes.SBE38temperature.value] = (
                    int(hex_segment[n : n + HEX_LENGTH["SBE38temperature"]], 16) / 100000 - 10
                )
                n += HEX_LENGTH["SBE38temperature"]

            if sensor == Sensors.WETLABS:
                results[HexDataTypes.wetlabs0.value] = int(
                    hex_segment[n : n + HEX_LENGTH["wetlabsSingleSensor"]], 16
                )
                n += HEX_LENGTH["wetlabsSingleSensor"]

                results[HexDataTypes.wetlabs1.value] = int(
                    hex_segment[n : n + HEX_LENGTH["wetlabsSingleSensor"]], 16
                )
                n += HEX_LENGTH["wetlabsSingleSensor"]

                results[HexDataTypes.wetlabs2.value] = int(
                    hex_segment[n : n + HEX_LENGTH["wetlabsSingleSensor"]], 16
                )
                n += HEX_LENGTH["wetlabsSingleSensor"]

            if sensor == Sensors.GTD:
                results[HexDataTypes.GTDpressure.value] = (
                    int(hex_segment[n : n + HEX_LENGTH["GTDpressure"]], 16) / 10000
                )
                n += HEX_LENGTH["GTDpressure"]
                results[HexDataTypes.GTDtemperature.value] = (
                    int(hex_segment[n : n + HEX_LENGTH["GTDtemperature"]], 16) / 10000 - 10
                )
                n += HEX_LENGTH["GTDtemperature"]

            if sensor == Sensors.DualGTD:
                results[HexDataTypes.GTDpressure.value] = (
                    int(hex_segment[n : n + HEX_LENGTH["GTDpressure"]], 16) / 10000
                )
                n += HEX_LENGTH["GTDpressure"]
                results[HexDataTypes.GTDtemperature.value] = (
                    int(hex_segment[n : n + HEX_LENGTH["GTDtemperature"]], 16) / 10000 - 10
                )
                n += HEX_LENGTH["GTDtemperature"]
                results[HexDataTypes.GTDpressure2.value] = (
                    int(hex_segment[n : n + HEX_LENGTH["GTDpressure"]], 16) / 10000
                )
                n += HEX_LENGTH["GTDpressure"]
                results[HexDataTypes.GTDtemperature2.value] = (
                    int(hex_segment[n : n + HEX_LENGTH["GTDtemperature"]], 16) / 10000 - 10
                )
                n += HEX_LENGTH["GTDtemperature"]

            if sensor == Sensors.OPTODE:
                results[HexDataTypes.optodeOxygen.value] = (
                    int(hex_segment[n : n + HEX_LENGTH["optodeOxygen"]], 16) / 10000 - 10
                )
                n += HEX_LENGTH["optodeOxygen"]

            if sensor == Sensors.SBE63:
                results[HexDataTypes.SBE63phase.value] = (
                    int(hex_segment[n : n + HEX_LENGTH["SBE63phase"]], 16) / 100000 - 10
                )
                n += HEX_LENGTH["SBE63phase"]
                results[HexDataTypes.SBE63temperature.value] = (
                    int(hex_segment[n : n + HEX_LENGTH["SBE63temperature"]], 16) / 1000000 - 1
                )
                n += HEX_LENGTH["SBE63temperature"]

            # Extract NMEA Sensors
            if sensor == Sensors.nmeaLatitude:
                lat = read_nmea_coordinates(hex_segment[n : n + HEX_LENGTH["nmeaLatitude"]])
                results[HexDataTypes.nmeaLatitude.value] = lat
                n += HEX_LENGTH["nmeaLatitude"]
            if sensor == Sensors.nmeaLongitude:
                lon = read_nmea_coordinates(hex_segment[n : n + HEX_LENGTH["nmeaLongitude"]])
                results[HexDataTypes.nmeaLongitude.value] = lon
                n += HEX_LENGTH["nmeaLongitude"]
            if sensor == Sensors.statusAndSign:
                signs = read_status_sign(hex_segment[n : n + HEX_LENGTH["statusAndSign"]])
                results[HexDataTypes.nmeaLatitude.value] *= signs[0]
                results[HexDataTypes.nmeaLongitude.value] *= signs[1]
                n += HEX_LENGTH["statusAndSign"]
            if sensor == Sensors.nmeaTime:
                seconds_since_2000 = read_nmea_time(hex_segment[n : n + HEX_LENGTH["nmeaTime"]])
                # naive, no timezone conversion
                results[HexDataTypes.nmeaTime.value] = datetime(1970, 1, 1) + timedelta(
                    seconds=seconds_since_2000 + SECONDS_BETWEEN_EPOCH_AND_2000
                )
                n += HEX_LENGTH["nmeaTime"]

    if moored_mode:
        seconds_since_2000 = int(hex_segment[n : n + HEX_LENGTH["time"]], 16)
        results[HexDataTypes.dateTime.value] = datetime.fromtimestamp(
            seconds_since_2000 + SECONDS_BETWEEN_EPOCH_AND_2000
        )
        n += HEX_LENGTH["time"]

    # Validate hex length. Ensure length matches what is expected based
    # on enabled sensors and moored mode.
    if n != len(hex_segment.strip()):
        raise RuntimeWarning(
            "Hex string length does not match expectation based on enabled sensors and moored mode"
        )

    return results


# TODO: change this to be snake_case for TKIT-75
def read_SBE37SM_format_0(  # pylint: disable=invalid-name
    hex_segment: str = "",
    enabled_sensors: List[Sensors] | None = None,
    hex=hex,  # pylint: disable=redefined-builtin
) -> Dict[str, Union[int, float, datetime]]:
    """Converts a 37 family data hex string into engineering units.

    :param hex_segment: one line from a hex data file
    :param enabled_sensors: array of Sensors that are enabled. For 37
        this is always temperature, conductivity, pressure. Defaults to
        False
    :param hex: Deprecated, use hex_segment

    :return: the 37 family sensor values in engineering units that were
        extracted from the input hex string
    """
    if hex is not builtins.hex:
        warnings.warn("hex is deprecated, use hex_segment", DeprecationWarning)

    results: Dict[str, Union[int, float, datetime]] = {}
    n = 0
    results[HexDataTypes.temperature.value] = int(hex_segment[n : HEX_LENGTH["temperature"]], 16)
    n += HEX_LENGTH["temperature"]

    results[HexDataTypes.conductivity.value] = (
        int(hex_segment[n : n + HEX_LENGTH["conductivity"]], 16) / 256
    )
    n += HEX_LENGTH["conductivity"]

    if enabled_sensors and Sensors.SBE63 in enabled_sensors:
        results[HexDataTypes.SBE63phase.value] = (
            int(hex_segment[n : n + HEX_LENGTH["SBE63phase"]], 16) / 100000 - 10
        )
        n += HEX_LENGTH["SBE63phase"]
        results[HexDataTypes.SBE63temperature.value] = (
            int(hex_segment[n : n + HEX_LENGTH["SBE63temperature"]], 16) / 1000000 - 1
        )
        n += HEX_LENGTH["SBE63temperature"]

    if enabled_sensors and Sensors.Pressure in enabled_sensors:
        results[HexDataTypes.pressure.value] = int(hex_segment[n : n + HEX_LENGTH["pressure"]], 16)
        n += HEX_LENGTH["pressure"]
        result = int(hex_segment[n : n + HEX_LENGTH["temperatureCompensation"]], 16)
        results[HexDataTypes.temperatureCompensation.value] = result
        n += HEX_LENGTH["temperatureCompensation"]

    seconds_since_2000 = int(hex_segment[n : n + HEX_LENGTH["time"]], 16)
    results[HexDataTypes.dateTime.value] = datetime(1970, 1, 1) + timedelta(
        seconds=seconds_since_2000 + SECONDS_BETWEEN_EPOCH_AND_2000
    )
    n += HEX_LENGTH["time"]

    # Validate hex length
    if n != len(hex_segment.strip()):
        raise ValueError("Hex string length does not match expectation based on enabled sensors")

    return results


def read_nmea_coordinates(hex_segment: str):
    """Converts a 3 byte NMEA hex string to latitude or longitude

    :param hex_segment: 3 byte hex string
    :raises RuntimeWarning: raised if the hex string is the wrong length
    :return: latitude or longitide coordinate
    """
    if len(hex_segment) != 6:
        raise RuntimeWarning(
            f"Unknown Coordinate Format. Received Hex of length {len(hex_segment)}. "
            f"Should have received Hex of length {HEX_LENGTH['nmeaLongitude']}"
        )
    byte0 = int(hex_segment[0:2], 16)
    byte1 = int(hex_segment[2:4], 16)
    byte2 = int(hex_segment[4:6], 16)
    coordinate = (byte0 * 65536 + byte1 * 256 + byte2) / 50000
    return coordinate


def read_status_sign(hex_segment: str):
    """Converts a hex byte to the signs for NMEA latitude and longitude

    :param hex_segment: 1 byte hex string
    :raises RuntimeWarning: raised if the hex string is the wrong length
    :raises RuntimeWarning: raised when the signs are converted
        incorrectly
    :return: a list of two integers (1 or -1)
    """
    if len(hex_segment) != 2:
        raise RuntimeWarning("Unknown Status Format")
    integer = int(hex_segment, 16)
    binary = format(integer, "0>8b")
    signs = []
    if binary[0] == "0":
        signs.append(1)
    elif binary[0] == "1":
        signs.append(-1)

    if binary[1] == "0":
        signs.append(1)
    elif binary[1] == "1":
        signs.append(-1)
    if len(signs) != 2:
        raise RuntimeWarning("An error occured while processing Coordinate Signs")
    return signs


def read_nmea_time(hex_segment: str):
    """Convert an 8 byte hex string to the number of seconds since 2000

    :param hex_segment: an 8 byte hex string
    :raises RuntimeWarning: raised if the hex string is the wrong length
    :return: _description_
    """
    if len(hex_segment) != 8:
        raise RuntimeWarning("Unknown Time Format")
    byte0 = hex_segment[0:2]
    byte1 = hex_segment[2:4]
    byte2 = hex_segment[4:6]
    byte3 = hex_segment[6:8]
    reformatted = int(byte3 + byte2 + byte1 + byte0, 16)
    return reformatted


def frequency_from_3_bytes(hex_segment: str) -> float:
    """Convert 3 bytes to a frequency value in Hz.

    :param hex_segment: 3-byte (6-character) hex string
    :return: frequency in Hz
    """
    byte0 = int(hex_segment[0:2], 16)
    byte1 = int(hex_segment[2:4], 16)
    byte2 = int(hex_segment[4:6], 16)
    # SBE11plus V2 user manual version A, page 39
    frequency = byte0 * 256 + byte1 + byte2 / 256
    return frequency


def voltages_from_3_bytes(hex_segment: str) -> tuple[float, float]:
    """Convert 3 bytes to two voltage channels.

    Each voltage channel is 12 bits; adjacent channels share a byte.
    Channels are suppressed in pairs, 3 bytes at a time.

    :param hex_segment: 3-byte (6-character) hex string
    :return: tuple of two voltages (voltageA, voltageB)
    """
    decimal_a = int(hex_segment[0 : HEX_LENGTH["SBE911Voltage"]], 16)
    voltage_a = 5 * (1 - decimal_a / 4095)

    decimal_b = int(hex_segment[HEX_LENGTH["SBE911Voltage"] : HEX_LENGTH["SBE911Voltage"] * 2], 16)
    voltage_b = 5 * (1 - decimal_b / 4095)

    return voltage_a, voltage_b


def reverse_hex_bytes(bytes_str: str) -> str:
    """
    Reverse the ASCII hex byte ordering.

    :param bytes_str: The ASCII hex bytes that need reordering.

    :return: The reordered ASCII hex bytes.
    """
    if not bytes_str:
        raise ValueError("Nothing to reverse")

    if len(bytes_str) % 2 != 0:
        raise ValueError("Hex string length is not even")

    if len(bytes_str) == 2:
        return bytes_str

    swapped = ""
    num_bytes = len(bytes_str) // 2

    for idx in range(num_bytes - 1, -1, -1):
        swapped += bytes_str[idx * 2 : idx * 2 + 2]

    return swapped
