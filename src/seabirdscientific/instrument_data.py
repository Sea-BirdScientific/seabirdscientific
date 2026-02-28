"""Functions for processing instrument data."""

# Native imports
from enum import Enum
from datetime import datetime, timedelta
from logging import getLogger
from typing import List, Dict, Union
from pathlib import Path
import warnings

# Third-party imports
import numpy as np
import pandas as pd
import xarray as xr

# Sea-Bird imports
from seabirdscientific.utils import WarnAllMembersMeta


logger = getLogger(__name__)

COUNTS_TO_VOLTS = 13107
SECONDS_BETWEEN_EPOCH_AND_2000 = 946684800


"""Possible data types in hex files"""
HEX_TYPE_TEMPERATURE = "temperature"
HEX_TYPE_SECONDARY_TEMPERATURE = "secondary temperature"
HEX_TYPE_CONDUCTIVITY = "conductivity"
HEX_TYPE_SECONDARY_CONDUCTIVITY = "secondary conductivity"
HEX_TYPE_PRESSURE = "pressure"
HEX_TYPE_DIGIQUARTZ_PRESSURE = "digiquartz pressure"
HEX_TYPE_TEMPERATURE_COMPENSATION = "temperature compensation"
HEX_TYPE_EXTVOLT0 = "volt 0"
HEX_TYPE_EXTVOLT1 = "volt 1"
HEX_TYPE_EXTVOLT2 = "volt 2"
HEX_TYPE_EXTVOLT3 = "volt 3"
HEX_TYPE_EXTVOLT4 = "volt 4"
HEX_TYPE_EXTVOLT5 = "volt 5"
HEX_TYPE_EXTVOLT6 = "volt 6"
HEX_TYPE_EXTVOLT7 = "volt 7"
HEX_TYPE_SBE38_TEMPERATURE = "SBE38 temperature"
HEX_TYPE_WETLABS0 = "wetlabs - channel 0"
HEX_TYPE_WETLABS1 = "wetlabs - channel 1"
HEX_TYPE_WETLABS2 = "wetlabs - channel 2"
HEX_TYPE_GTD_PRESSURE = "GTD pressure"
HEX_TYPE_GTD_TEMPERATURE = "GTD temperature"
HEX_TYPE_GTD_PRESSURE2 = "GTD pressure - sensor 2"
HEX_TYPE_GTD_TEMPERATURE2 = "GTD temperature - sensor 2"
HEX_TYPE_OPTODE_OXYGEN = "optode oxygen"
HEX_TYPE_SBE63_PHASE = "SBE63 oxygen phase"
HEX_TYPE_SBE63_TEMPERATURE = "SBE63 oxygen temperature"
HEX_TYPE_DATE_TIME = "date time"
HEX_TYPE_NMEA_TIME = "NMEA Date Time"
HEX_TYPE_NMEA_LATITUDE = "NMEA Latitude"
HEX_TYPE_NMEA_LONGITUDE = "NMEA Longitude"
HEX_TYPE_NMEA_DEPTH = "nmea depth"
HEX_TYPE_STATUS_SIGN = "status and sign"
HEX_TYPE_VRS_INTERNAL = "vrs internal"
HEX_TYPE_VRS_EXTERNAL = "vrs external"
HEX_TYPE_PH_TEMPERATURE = "ph temperature"
HEX_TYPE_VK = "vk"
HEX_TYPE_IB = "ib"
HEX_TYPE_IK = "ik"
HEX_TYPE_RELATIVE_HUMIDITY = "relative humidity"
HEX_TYPE_INTERNAL_TEMPERATURE = "internal temperature"
HEX_TYPE_ERROR_FLAG = "error flag"
HEX_TYPE_SURFACE_PAR = "surface par"
HEX_TYPE_SBE911_PUMP_STATUS = "SBE911 pump status"
HEX_TYPE_SBE911_BOTTOM_CONTACT_STATUS = "SBE911 bottom contact status"
HEX_TYPE_SBE911_CONFIRM_STATUS = "SBE911 confirm status"
HEX_TYPE_SBE911_MODEM_STATUS = "SBE911 modem status"
HEX_TYPE_DATA_INTEGRITY = "data integrity"
HEX_TYPE_SYSTEM_TIME = "system time"


"""Possible lengths for hex data types"""
HEX_LEN_TEMPERATURE = 6
HEX_LEN_CONDUCTIVITY = 6
HEX_LEN_PRESSURE = 6
HEX_LEN_TEMPERATURE_COMPENSATION = 4
HEX_LEN_VOLTAGE = 4
HEX_LEN_SBE911_VOLTAGE = 3
HEX_LEN_SBE911_SURFACE_PAR = 3
HEX_LEN_SBE911_TEMPERATURE_COMPENSATION = 3
HEX_LEN_SBE911_STATUS = 1
HEX_LEN_SBE911_DATA_INTEGRITY = 2
HEX_LEN_WETLABS_SINGLE_SENSOR = 4
HEX_LEN_GTD_PRESSURE = 8
HEX_LEN_OPTODE_OXYGEN = 6
HEX_LEN_SBE63_PHASE = 6
HEX_LEN_SEAOWL_CHANNEL = 4
HEX_LEN_DATE_TIME = 8
HEX_LEN_VRS_EXTERNAL = 6
HEX_LEN_VRS_INTERNAL = 6
HEX_LEN_VK = 6
HEX_LEN_IB = 6
HEX_LEN_IK = 6
HEX_LEN_RELATIVE_HUMIDITY = 3
HEX_LEN_INTERNAL_TEMPERATURE = 3
HEX_LEN_ERROR_FLAG = 4
HEX_LEN_NMEA_LATITUDE = 6
HEX_LEN_NMEA_LONGITUDE = 6
HEX_LEN_NMEA_TIME = 8
HEX_LEN_NMEA_LOCATION = 14
HEX_LEN_NMEA_STATUS_AND_SIGN = 2
HEX_LEN_NMEA_DEPTH = 6
HEX_LEN_SYSTEM_TIME = 8


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
    SBE19Plus = "19plus"
    SBE16Plus = "16plus"
    SBE39Plus = "39plus"
    SBE39PlusIM = "39plus-IM"
    SBE911Plus = "911plus"
    SeaFET2 = "SeaFET2"
    SeapHox2 = "SeapHox2"
    HydroCAT = "HydroCAT"
    HydroCATODO = "HydroCAT-ODO"


class Sensors(Enum):
    """Available sensors to read hex data from"""

    Temperature = "Temperature"
    SecondaryTemperature = "SecondaryTemperature"
    Conductivity = "Conductivity"
    SecondaryConductivity = "SecondaryConductivity"
    Pressure = "Pressure"
    ExtVolt0 = "ExtVolt0"
    ExtVolt1 = "ExtVolt1"
    ExtVolt2 = "ExtVolt2"
    ExtVolt3 = "ExtVolt3"
    ExtVolt4 = "ExtVolt4"
    ExtVolt5 = "ExtVolt5"
    ExtVolt6 = "ExtVolt6"
    ExtVolt7 = "ExtVolt7"
    WETLABS = "WETLABS"
    GTD = "GTD"
    DualGTD = "DualGTD"
    OPTODE = "OPTODE"
    SBE63 = "SBE63"
    SBE38 = "SBE38"
    SeaFET = "SeaFET"
    SPAR = "SPAR"
    nmeaLatitude = "nmeaLatitude"
    nmeaLongitude = "nmeaLongitude"
    statusAndSign = "StatusAndSign"
    nmeaTime = "nmeaTime"
    nmeaLocation = "nmeaLocation"
    nmeaDepth = "nmeaDepth"
    SystemTime = "systemTime"


class HexDataTypes(Enum, metaclass=WarnAllMembersMeta):
    """Possible data types in hex files.
    Deprecated. Use HEX_TYPE_* constants."""

    temperature = "temperature"
    secondaryTemperature = "secondary temperature"
    conductivity = "conductivity"
    secondaryConductivity = "secondary conductivity"
    pressure = "pressure"
    digiquartzPressure = "digiquartz pressure"
    temperatureCompensation = "temperature compensation"
    ExtVolt0 = "volt 0"
    ExtVolt1 = "volt 1"
    ExtVolt2 = "volt 2"
    ExtVolt3 = "volt 3"
    ExtVolt4 = "volt 4"
    ExtVolt5 = "volt 5"
    ExtVolt6 = "volt 6"
    ExtVolt7 = "volt 7"
    SBE38temperature = "SBE38 temperature"
    wetlabs0 = "wetlabs - channel 0"
    wetlabs1 = "wetlabs - channel 1"
    wetlabs2 = "wetlabs - channel 2"
    GTDpressure = "GTD pressure"
    GTDtemperature = "GTD temperature"
    GTDpressure2 = "GTD pressure - sensor 2"
    GTDtemperature2 = "GTD temperature - sensor 2"
    optodeOxygen = "optode oxygen"
    SBE63phase = "SBE63 oxygen phase"
    SBE63temperature = "SBE63 oxygen temperature"
    dateTime = "date time"
    nmeaTime = "NMEA Date Time"
    nmeaLatitude = "NMEA Latitude"
    nmeaLongitude = "NMEA Longitude"
    nmeaDepth = "nmea depth"
    statusAndSign = "status and sign"
    vrsInternal = "vrs internal"
    vrsExternal = "vrs external"
    pHtemperature = "ph temperature"
    vk = "vk"
    ib = "ib"
    ik = "ik"
    relativeHumidity = "relative humidity"
    internalTemperature = "internal temperature"
    errorFlag = "error flag"
    surfacePAR = "surface par"
    SBE911PumpStatus = "SBE911 pump status"
    SBE911BottomContactStatus = "SBE911 bottom contact status"
    SBE911ConfirmStatus = "SBE911 confirm status"
    SBE911ModemStatus = "SBE911 modem status"
    dataIntegrity = "data integrity"
    systemTime = "system time"


def read_cnv_file(filepath: Union[Path, str]) -> xr.Dataset:
    """Import the data from a .cnv file and put it into an xarray
    Dataset. Duplicate varioable names will have a number appended. For
    example, the second "depSM" becomes "depSM_1".

    :param filepath: the path to the .cnv file to be imported

    :return: the imported data from the .cnv file

    """

    dataset = xr.Dataset({}, attrs={"file_name": Path(filepath).name})

    total_samples = 0
    data_lines = []

    logger.info("Unpacking instrument data from file: %s", filepath)

    with open(filepath, mode="r") as cnv:
        for n_line, line in enumerate(cnv):
            if line.startswith("# nvalues = "):
                total_samples = int(line[line.find("= ") + 2 : line.find("\n")])

            elif line.startswith("# name "):
                name = line[line.find("= ") + 2 : line.find(":")]
                left_bracket = line.find(" [")
                if left_bracket > 0:
                    long_name = line[line.find(": ") + 2 : left_bracket]
                    units = line[line.find("[") + 1 : line.find("]")]
                else:
                    long_name = line[line.find(": ") + 2 : line.find("\n")].strip()
                    units = ""

                safe_name = name
                # check for label collision
                # if collision, increment until find available key
                n_name = 1
                while safe_name in list(dataset.data_vars):
                    safe_name = f"{name}_{n_name}"
                    n_name += 1

                if safe_name != name:
                    logger.warning(f'Duplicate measurand "{name}" will use key "{safe_name}"')

                data_array = xr.DataArray(
                    data=np.zeros(total_samples),
                    dims=["sample"],
                    coords={"sample": np.arange(total_samples)},
                    attrs={
                        "sbs_name": name,
                        "long_name": long_name,
                        "units": units,
                    },
                )

                dataset[safe_name] = data_array

            elif line.startswith("# interval = "):
                interval = float(
                    line[line.find(": ") + 2 : line.find("\n")]
                )  # TODO: fix for minutes, hours, etc
                dataset.attrs["sample_interval"] = interval

            elif line.startswith("# start_time = "):
                end = min(line.find("\n"), line.find(" ["))
                date_string = line[line.find("= ") + 2 : end]
                start_time = datetime.strptime(date_string, "%b %d %Y %H:%M:%S")
                dataset.attrs["start_time"] = start_time

            elif line.startswith("*END*"):
                data_lines = cnv.readlines()
                break

    np_data = np.array(
        [np.fromstring(dl.replace("-", " -").replace("e -", "e-"), sep=" ") for dl in data_lines]
    )

    for n, measurand in enumerate(list(dataset.data_vars)):
        dataset[measurand].data = np_data[:, n]

    return dataset


def read_hex_file(
    filepath: Union[Path, str],
    instrument_type: InstrumentType,
    enabled_sensors: List[Sensors] = [],
    moored_mode=False,
    is_shallow=True,
    frequency_channels_suppressed=0,
    voltage_words_suppressed=0,
) -> xr.Dataset:
    """Reads a .hex file

    :param filepath: path to the .hex file
    :param instrument_type: the instrument that generated the .hex file
    :param enabled_sensors: list of sensors that were enabled on the
        instrument
    :param moored_mode: whether the instrument was in moored or profiling
        mode, defaults to False
    :param is_shallow: boolean for deep or shallow seafet,
    :param frequency_channels_suppressed: number of SBE911 requency
        channels supressed,
    :param voltage_words_suppressed: number of SBE911 voltage channels
        suppressed,
    :return: an xarray Dataset with the hex data
    """

    data_lines = []

    with open(filepath, mode="r") as file:
        for line in file:
            if line.startswith("*END*"):
                data_lines = file.readlines()

    dataset = xr.Dataset({}, attrs={"file_name": Path(filepath).name})

    hex_data = read_hex(
        instrument_type,
        data_lines[0],
        enabled_sensors,
        moored_mode,
        is_shallow,
        frequency_channels_suppressed,
        voltage_words_suppressed,
    )
    keys = hex_data.keys()
    dataset.update(_preallocate_dataset(hex_data, len(data_lines)))
    np_data = np.empty((len(data_lines), len(keys)), dtype=object)

    for n, line in enumerate(data_lines):
        hex_data = read_hex(
            instrument_type,
            line,
            enabled_sensors,
            moored_mode,
            is_shallow,
            frequency_channels_suppressed,
            voltage_words_suppressed,
        )

        np_data[n] = [hex_data[k] for k in keys]

    for n, measurand in enumerate(keys):
        dtype = np.dtype(type(hex_data[measurand]))
        dataset[measurand].data = np_data[:, n].astype(dtype)

    return dataset


def _preallocate_dataset(
    hex_data: dict,
    data_length: int,
) -> xr.Dataset:
    """Creates an xarray dataset prefilled with zeros or arbitrary dates
    to replace with instrument data

    :param hex_data: a sample from a hex file converted to a dict
    :param data_length: the number of samples in the hex file

    :return: a dataset fill of zeros or arbitrary dates
    """
    dataset = xr.Dataset()
    for key, value in hex_data.items():
        if isinstance(value, datetime):
            # fill with arbitrary dates converted to regular datetimes
            date_range = pd.date_range(start="2000-01-01", end="2000-01-02", periods=data_length)
            empty_data = date_range.to_pydatetime().tolist()
        else:
            empty_data = np.zeros(data_length)

        data_array = xr.DataArray(
            data=empty_data,
            dims=["sample"],
            coords={"sample": np.arange(len(empty_data))},
            attrs={
                "sbs_name": key,
                "long_name": "",
                "units": "",
            },
        )
        dataset[key] = data_array

    return dataset


def read_hex(
    instrument_type: InstrumentType,
    hex_segment: str = "",
    enabled_sensors: Union[List[Sensors], None] = None,
    moored_mode=False,
    is_shallow=True,
    frequency_channels_suppressed=0,
    voltage_words_suppressed=0,
) -> dict:
    """Converts an instrument data hex string into engineering units.

    :param instrument_type: determines how units are converted
    :param hex_segment: one line from a hex data file
    :param enabled_sensors: mooredMode parses time for 19plus in moored
        mode if true
    :param moored_mode: array of Sensors that are enabled. For 37 this
        is always temperature, conductivity, pressure. Defaults to False

    :return: the sensor values in engineering units that were extracted
        from the input hex string
    """

    if instrument_type == InstrumentType.SBE19Plus:
        return read_sbe19plus_format_0(hex_segment, enabled_sensors, moored_mode)

    elif instrument_type == InstrumentType.SBE16Plus:
        return read_sbe19plus_format_0(hex_segment, enabled_sensors, True)

    elif instrument_type in [
        InstrumentType.SBE37IM,
        InstrumentType.SBE37IMP,
        InstrumentType.SBE37IMPODO,
        InstrumentType.SBE37IM,
        InstrumentType.SBE37SI,
        InstrumentType.SBE37SIP,
        InstrumentType.SBE37SM,
        InstrumentType.SBE37SMP,
        InstrumentType.SBE37SMPODO,
        InstrumentType.HydroCAT,
        InstrumentType.HydroCATODO,
    ]:
        return read_sbe37sm_format_0(hex_segment, enabled_sensors)

    elif instrument_type in [InstrumentType.SBE39Plus, InstrumentType.SBE39PlusIM]:
        return read_sbe39plus_format_0(hex_segment, enabled_sensors)

    elif instrument_type in [InstrumentType.SeaFET2, InstrumentType.SeapHox2]:
        return read_seafet_format_0(hex_segment, instrument_type, is_shallow)

    elif instrument_type == InstrumentType.SBE911Plus:
        return read_sbe911plus_format_0(
            hex_segment,
            enabled_sensors,
            frequency_channels_suppressed,
            voltage_words_suppressed,
        )

    else:
        raise ValueError(f"Unsupported instrument type: {instrument_type}")

    return {}


def read_sbe39plus_format_0(
    hex_segment: str = "",
    enabled_sensors: Union[List[Sensors], None] = None,
) -> Dict[str, Union[int, float, datetime]]:
    """Converts a 39plus data hex string into engineering units.

    :param hex_segment: one line from a hex data file
    :param enabled_sensors: array of Sensors that are enabled

    :return: the 39plus sensor values in engineering units that were
        extracted from the input hex string
    """

    results: Dict[str, Union[int, float, datetime]] = {}
    n = 0

    # Datetime (naive, no timezone conversion)
    seconds_since_2000 = int(hex_segment[n : n + HEX_LEN_DATE_TIME], 16)
    results[HEX_TYPE_DATE_TIME] = datetime(1970, 1, 1) + timedelta(
        seconds=seconds_since_2000 + SECONDS_BETWEEN_EPOCH_AND_2000
    )
    n += HEX_LEN_DATE_TIME

    # Temperature
    results[HEX_TYPE_TEMPERATURE] = int(hex_segment[n : n + HEX_LEN_TEMPERATURE], 16)
    n += HEX_LEN_TEMPERATURE

    # Pressure and temperature compensation
    if enabled_sensors and Sensors.Pressure in enabled_sensors:
        results[HEX_TYPE_PRESSURE] = int(hex_segment[n : n + HEX_LEN_PRESSURE], 16)
        n += HEX_LEN_PRESSURE

        temp_comp = int(hex_segment[n : n + HEX_LEN_TEMPERATURE_COMPENSATION], 16)
        results[HEX_TYPE_TEMPERATURE_COMPENSATION] = temp_comp
        n += HEX_LEN_TEMPERATURE_COMPENSATION

    # Validate hex length
    if n != len(hex_segment.strip()):
        raise ValueError("Hex string length does not match expectation based on enabled sensors")

    return results


def read_SBE39plus_format_0(*args, **kwargs):
    warnings.warn("Deprecated, use read_sbe39plus_format_0", DeprecationWarning)
    return read_sbe39plus_format_0(*args, **kwargs)


def read_seafet_format_0(
    hex_segment: str,
    instrument_type: InstrumentType,
    is_shallow: bool = True,
) -> Dict[str, Union[int, float, datetime]]:
    """Converts a SeaFET2 or SeapHox2 hex string into engineering units.

    :param hex_segment: one line from a hex data file
    :param instrument_type: InstrumentType.SeaFET2 or InstrumentType.SeapHox2
    :param is_shallow: if True, include internal pH and pH reference temperature

    :return: sensor values in engineering units extracted from the hex string
    """

    if instrument_type not in (InstrumentType.SeaFET2, InstrumentType.SeapHox2):
        raise ValueError(
            f"In read_seafet_format_0 {instrument_type} is not recognized as a SeaFET2 or SeapHox2"
        )

    results: Dict[str, Union[int, float, datetime]] = {}
    n = 0

    # SeapHox2 specific values
    if instrument_type == InstrumentType.SeapHox2:
        results[HEX_TYPE_TEMPERATURE] = int(hex_segment[n : n + HEX_LEN_TEMPERATURE], 16)
        n += HEX_LEN_TEMPERATURE

        results[HEX_TYPE_CONDUCTIVITY] = int(hex_segment[n : n + HEX_LEN_CONDUCTIVITY], 16) / 256
        n += HEX_LEN_CONDUCTIVITY

        results[HEX_TYPE_PRESSURE] = int(hex_segment[n : n + HEX_LEN_PRESSURE], 16)
        n += HEX_LEN_PRESSURE

        results[HEX_TYPE_TEMPERATURE_COMPENSATION] = int(
            hex_segment[n : n + HEX_LEN_TEMPERATURE_COMPENSATION], 16
        )
        n += HEX_LEN_TEMPERATURE_COMPENSATION

        results[HEX_TYPE_SBE63_PHASE] = (
            int(hex_segment[n : n + HEX_LEN_SBE63_PHASE], 16) / 100000 - 10
        )
        n += HEX_LEN_SBE63_PHASE

        results[HEX_TYPE_SBE63_TEMPERATURE] = (
            int(hex_segment[n : n + HEX_LEN_TEMPERATURE], 16) / 1000000 - 1
        )
        n += HEX_LEN_TEMPERATURE

    # External pH
    results[HEX_TYPE_VRS_EXTERNAL] = int(hex_segment[n : n + HEX_LEN_VRS_EXTERNAL], 16)
    n += HEX_LEN_VRS_EXTERNAL

    if is_shallow:
        # Internal pH
        results[HEX_TYPE_VRS_INTERNAL] = int(hex_segment[n : n + HEX_LEN_VRS_INTERNAL], 16)
        n += HEX_LEN_VRS_INTERNAL

        # pH reference temperature
        results[HEX_TYPE_PH_TEMPERATURE] = int(hex_segment[n : n + HEX_LEN_TEMPERATURE], 16)
        n += HEX_LEN_TEMPERATURE

    # Other sensors
    results[HEX_TYPE_VK] = int(hex_segment[n : n + HEX_LEN_VK], 16)
    n += HEX_LEN_VK

    results[HEX_TYPE_IB] = int(hex_segment[n : n + HEX_LEN_IB], 16)
    n += HEX_LEN_IB

    results[HEX_TYPE_IK] = int(hex_segment[n : n + HEX_LEN_IK], 16)
    n += HEX_LEN_IK

    results[HEX_TYPE_RELATIVE_HUMIDITY] = int(
        hex_segment[n : n + HEX_LEN_RELATIVE_HUMIDITY] + "0", 16
    )
    n += HEX_LEN_RELATIVE_HUMIDITY

    results[HEX_TYPE_INTERNAL_TEMPERATURE] = int(
        hex_segment[n : n + HEX_LEN_INTERNAL_TEMPERATURE] + "0", 16
    )
    n += HEX_LEN_INTERNAL_TEMPERATURE

    # Datetime (naive, no timezone conversion)
    seconds_since_2000 = int(hex_segment[n : n + HEX_LEN_DATE_TIME], 16)
    results[HEX_TYPE_DATE_TIME] = datetime(2000, 1, 1) + timedelta(seconds=seconds_since_2000)
    n += HEX_LEN_DATE_TIME

    # Error flag
    results[HEX_TYPE_ERROR_FLAG] = int(hex_segment[n : n + HEX_LEN_ERROR_FLAG], 16)
    n += HEX_LEN_ERROR_FLAG

    # Validate hex length
    if n != len(hex_segment.strip()):
        raise ValueError(
            "Hex string length does not match expectation based on instrument_type and is_shallow"
        )

    return results


def read_sbe911plus_format_0(
    hex_segment: str = "",
    enabled_sensors: Union[List[Sensors], None] = None,
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
    results[HEX_TYPE_TEMPERATURE] = frequency_from_3_bytes(
        hex_segment[n : n + HEX_LEN_TEMPERATURE]
    )
    n += HEX_LEN_TEMPERATURE

    # Conductivity
    results[HEX_TYPE_CONDUCTIVITY] = frequency_from_3_bytes(
        hex_segment[n : n + HEX_LEN_CONDUCTIVITY]
    )
    n += HEX_LEN_CONDUCTIVITY

    # Digiquartz Pressure
    results[HEX_TYPE_DIGIQUARTZ_PRESSURE] = frequency_from_3_bytes(
        hex_segment[n : n + HEX_LEN_PRESSURE]
    )
    n += HEX_LEN_PRESSURE

    # Secondary temperature
    if Sensors.SecondaryTemperature in enabled_sensors:
        results[HEX_TYPE_SECONDARY_TEMPERATURE] = frequency_from_3_bytes(
            hex_segment[n : n + HEX_LEN_TEMPERATURE]
        )
        n += HEX_LEN_TEMPERATURE
    elif frequency_channels_suppressed <= 1:
        n += HEX_LEN_TEMPERATURE

    # Secondary conductivity
    if Sensors.SecondaryConductivity in enabled_sensors:
        results[HEX_TYPE_SECONDARY_CONDUCTIVITY] = frequency_from_3_bytes(
            hex_segment[n : n + HEX_LEN_CONDUCTIVITY]
        )
        n += HEX_LEN_CONDUCTIVITY
    elif frequency_channels_suppressed == 0:
        n += HEX_LEN_CONDUCTIVITY

    # Voltage channels (suppressed in pairs)
    if Sensors.ExtVolt0 in enabled_sensors or Sensors.ExtVolt1 in enabled_sensors:
        results[HEX_TYPE_EXTVOLT0], results[HEX_TYPE_EXTVOLT1] = voltages_from_3_bytes(
            hex_segment[n : n + HEX_LEN_SBE911_VOLTAGE * 2]
        )
        n += HEX_LEN_SBE911_VOLTAGE * 2
    elif voltage_words_suppressed <= 3:
        # NotInUse volt channel that was not suppressed
        n += HEX_LEN_SBE911_VOLTAGE * 2

    if Sensors.ExtVolt2 in enabled_sensors or Sensors.ExtVolt3 in enabled_sensors:
        results[HEX_TYPE_EXTVOLT2], results[HEX_TYPE_EXTVOLT3] = voltages_from_3_bytes(
            hex_segment[n : n + HEX_LEN_SBE911_VOLTAGE * 2]
        )
        n += HEX_LEN_SBE911_VOLTAGE * 2
    elif voltage_words_suppressed <= 2:
        # NotInUse volt channel that was not suppressed
        n += HEX_LEN_SBE911_VOLTAGE * 2

    if Sensors.ExtVolt4 in enabled_sensors or Sensors.ExtVolt5 in enabled_sensors:
        results[HEX_TYPE_EXTVOLT4], results[HEX_TYPE_EXTVOLT5] = voltages_from_3_bytes(
            hex_segment[n : n + HEX_LEN_SBE911_VOLTAGE * 2]
        )
        n += HEX_LEN_SBE911_VOLTAGE * 2
    elif voltage_words_suppressed <= 1:
        # NotInUse volt channel that was not suppressed
        n += HEX_LEN_SBE911_VOLTAGE * 2

    if Sensors.ExtVolt6 in enabled_sensors or Sensors.ExtVolt7 in enabled_sensors:
        results[HEX_TYPE_EXTVOLT6], results[HEX_TYPE_EXTVOLT7] = voltages_from_3_bytes(
            hex_segment[n : n + HEX_LEN_SBE911_VOLTAGE * 2]
        )
        n += HEX_LEN_SBE911_VOLTAGE * 2
    elif voltage_words_suppressed == 0:
        # NotInUse volt channel that was not suppressed
        n += HEX_LEN_SBE911_VOLTAGE * 2

    # Surface PAR
    if Sensors.SPAR in enabled_sensors:
        n += 3  # unused bits
        results[HEX_TYPE_SURFACE_PAR] = (
            int(hex_segment[n : n + HEX_LEN_SBE911_SURFACE_PAR], 16) / 819
        )
        n += HEX_LEN_SBE911_SURFACE_PAR

    # NMEA Location
    if Sensors.nmeaLocation in enabled_sensors:
        hex_loc = hex_segment[n : n + HEX_LEN_NMEA_LOCATION]

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

        (
            results[HEX_TYPE_NMEA_LATITUDE],
            results[HEX_TYPE_NMEA_LONGITUDE],
        ) = get_lat_lon(hex_loc)
        n += HEX_LEN_NMEA_LOCATION

    # NMEA Depth
    if Sensors.nmeaDepth in enabled_sensors:
        results[HEX_TYPE_NMEA_DEPTH] = int(hex_segment[n : n + HEX_LEN_NMEA_DEPTH], 16)
        n += HEX_LEN_NMEA_DEPTH

    # NMEA Time
    if Sensors.nmeaTime in enabled_sensors:
        seconds_since_2000 = int(reverse_hex_bytes(hex_segment[n : n + HEX_LEN_NMEA_TIME]), 16)
        results[HEX_TYPE_NMEA_TIME] = datetime(2000, 1, 1) + timedelta(seconds=seconds_since_2000)
        n += HEX_LEN_NMEA_TIME

    # Temperature compensation
    results[HEX_TYPE_TEMPERATURE_COMPENSATION] = int(
        hex_segment[n : n + HEX_LEN_SBE911_TEMPERATURE_COMPENSATION], 16
    )
    n += HEX_LEN_SBE911_TEMPERATURE_COMPENSATION

    # Status bits
    status_bin = format(int(hex_segment[n : n + HEX_LEN_SBE911_STATUS], 16), "04b")
    results[HEX_TYPE_SBE911_PUMP_STATUS] = int(status_bin[0])
    results[HEX_TYPE_SBE911_BOTTOM_CONTACT_STATUS] = int(status_bin[1])
    results[HEX_TYPE_SBE911_CONFIRM_STATUS] = int(status_bin[2])
    results[HEX_TYPE_SBE911_MODEM_STATUS] = int(status_bin[3])
    n += HEX_LEN_SBE911_STATUS

    # Data integrity
    results[HEX_TYPE_DATA_INTEGRITY] = int(hex_segment[n : n + HEX_LEN_SBE911_DATA_INTEGRITY], 16)
    n += HEX_LEN_SBE911_DATA_INTEGRITY

    # System time
    if Sensors.SystemTime in enabled_sensors:
        seconds_since_1970 = int(reverse_hex_bytes(hex_segment[n : n + HEX_LEN_SYSTEM_TIME]), 16)
        results[HEX_TYPE_SYSTEM_TIME] = datetime(1970, 1, 1) + timedelta(
            seconds=seconds_since_1970
        )
        n += HEX_LEN_SYSTEM_TIME

    # Validate hex length
    if n != len(hex_segment.strip()):
        raise ValueError("Hex string length does not match expectation based on enabled sensors")

    return results


def read_SBE911plus_format_0(*args, **kwargs):
    warnings.warn("Deprecated, use read_sbe911plus_format_0", DeprecationWarning)
    return read_sbe911plus_format_0(*args, **kwargs)


def read_sbe19plus_format_0(
    hex_segment: str = "",
    enabled_sensors: Union[List[Sensors], None] = None,
    moored_mode=False,
) -> Dict[str, Union[float, datetime]]:
    """Converts a 19plus V2 data hex string into engineering units.

    :param hex_segment: one line from a hex data file
    :param enabled_sensors: array of Sensors that are enabled. For 37
        this is always temperature, conductivity, pressure. Defaults to
        None
    :param moored_mode: parses time for 19plus in moored mode if true.
        Defautls to False

    :return: the 19plus V2 sensor values in engineering units that were
            extracted from the input hex string

    :raises RuntimeWarning: if the hex string length does not match the
        expected length
    """

    results: Dict[str, Union[int, float, datetime]] = {}
    n = 0
    for sensor in Sensors:
        if enabled_sensors and sensor in enabled_sensors:
            if sensor == Sensors.Temperature:
                results[HEX_TYPE_TEMPERATURE] = int(hex_segment[n:HEX_LEN_TEMPERATURE], 16)
                n += HEX_LEN_TEMPERATURE

            if sensor == Sensors.Conductivity:
                results[HEX_TYPE_CONDUCTIVITY] = (
                    int(hex_segment[n : n + HEX_LEN_CONDUCTIVITY], 16) / 256
                )
                n += HEX_LEN_CONDUCTIVITY

            if sensor == Sensors.Pressure:  # TODO: add conversion for quartz pressure sensors
                results[HEX_TYPE_PRESSURE] = int(hex_segment[n : n + HEX_LEN_PRESSURE], 16)
                n += HEX_LEN_PRESSURE
                result = (
                    int(hex_segment[n : n + HEX_LEN_TEMPERATURE_COMPENSATION], 16)
                    / COUNTS_TO_VOLTS
                )
                results[HEX_TYPE_TEMPERATURE_COMPENSATION] = result
                n += HEX_LEN_TEMPERATURE_COMPENSATION

            if sensor == Sensors.ExtVolt0:
                results[HEX_TYPE_EXTVOLT0] = (
                    int(hex_segment[n : n + HEX_LEN_VOLTAGE], 16) / COUNTS_TO_VOLTS
                )
                n += HEX_LEN_VOLTAGE

            if sensor == Sensors.ExtVolt1:
                results[HEX_TYPE_EXTVOLT1] = (
                    int(hex_segment[n : n + HEX_LEN_VOLTAGE], 16) / COUNTS_TO_VOLTS
                )
                n += HEX_LEN_VOLTAGE

            if sensor == Sensors.ExtVolt2:
                results[HEX_TYPE_EXTVOLT2] = (
                    int(hex_segment[n : n + HEX_LEN_VOLTAGE], 16) / COUNTS_TO_VOLTS
                )
                n += HEX_LEN_VOLTAGE

            if sensor == Sensors.ExtVolt3:
                results[HEX_TYPE_EXTVOLT3] = (
                    int(hex_segment[n : n + HEX_LEN_VOLTAGE], 16) / COUNTS_TO_VOLTS
                )
                n += HEX_LEN_VOLTAGE

            if sensor == Sensors.ExtVolt4:
                results[HEX_TYPE_EXTVOLT4] = (
                    int(hex_segment[n : n + HEX_LEN_VOLTAGE], 16) / COUNTS_TO_VOLTS
                )
                n += HEX_LEN_VOLTAGE

            if sensor == Sensors.ExtVolt5:
                results[HEX_TYPE_EXTVOLT5] = (
                    int(hex_segment[n : n + HEX_LEN_VOLTAGE], 16) / COUNTS_TO_VOLTS
                )
                n += HEX_LEN_VOLTAGE

            if sensor == Sensors.SBE38:
                results[HEX_TYPE_SBE38_TEMPERATURE] = (
                    int(hex_segment[n : n + HEX_LEN_TEMPERATURE], 16) / 100000 - 10
                )
                n += HEX_LEN_TEMPERATURE

            if sensor == Sensors.WETLABS:
                results[HEX_TYPE_WETLABS0] = int(
                    hex_segment[n : n + HEX_LEN_WETLABS_SINGLE_SENSOR], 16
                )
                n += HEX_LEN_WETLABS_SINGLE_SENSOR

                results[HEX_TYPE_WETLABS1] = int(
                    hex_segment[n : n + HEX_LEN_WETLABS_SINGLE_SENSOR], 16
                )
                n += HEX_LEN_WETLABS_SINGLE_SENSOR

                results[HEX_TYPE_WETLABS2] = int(
                    hex_segment[n : n + HEX_LEN_WETLABS_SINGLE_SENSOR], 16
                )
                n += HEX_LEN_WETLABS_SINGLE_SENSOR

            if sensor == Sensors.GTD:
                results[HEX_TYPE_GTD_PRESSURE] = (
                    int(hex_segment[n : n + HEX_LEN_GTD_PRESSURE], 16) / 10000
                )
                n += HEX_LEN_GTD_PRESSURE

                results[HEX_TYPE_GTD_TEMPERATURE] = (
                    int(hex_segment[n : n + HEX_LEN_TEMPERATURE], 16) / 10000 - 10
                )
                n += HEX_LEN_TEMPERATURE

            if sensor == Sensors.DualGTD:
                results[HEX_TYPE_GTD_PRESSURE] = (
                    int(hex_segment[n : n + HEX_LEN_GTD_PRESSURE], 16) / 10000
                )
                n += HEX_LEN_GTD_PRESSURE

                results[HEX_TYPE_GTD_TEMPERATURE] = (
                    int(hex_segment[n : n + HEX_LEN_TEMPERATURE], 16) / 10000 - 10
                )
                n += HEX_LEN_TEMPERATURE

                results[HEX_TYPE_GTD_PRESSURE2] = (
                    int(hex_segment[n : n + HEX_LEN_GTD_PRESSURE], 16) / 10000
                )
                n += HEX_LEN_GTD_PRESSURE

                results[HEX_TYPE_GTD_TEMPERATURE2] = (
                    int(hex_segment[n : n + HEX_LEN_TEMPERATURE], 16) / 10000 - 10
                )
                n += HEX_LEN_TEMPERATURE

            if sensor == Sensors.OPTODE:
                results[HEX_TYPE_OPTODE_OXYGEN] = (
                    int(hex_segment[n : n + HEX_LEN_OPTODE_OXYGEN], 16) / 10000 - 10
                )
                n += HEX_LEN_OPTODE_OXYGEN

            if sensor == Sensors.SBE63:
                results[HEX_TYPE_SBE63_PHASE] = (
                    int(hex_segment[n : n + HEX_LEN_SBE63_PHASE], 16) / 100000 - 10
                )
                n += HEX_LEN_SBE63_PHASE
                results[HEX_TYPE_SBE63_TEMPERATURE] = (
                    int(hex_segment[n : n + HEX_LEN_TEMPERATURE], 16) / 1000000 - 1
                )
                n += HEX_LEN_TEMPERATURE

            # Extract NMEA Sensors
            if sensor == Sensors.nmeaLatitude:
                lat = read_nmea_coordinates(hex_segment[n : n + HEX_LEN_NMEA_LATITUDE])
                results[HEX_TYPE_NMEA_LATITUDE] = lat
                n += HEX_LEN_NMEA_LATITUDE

            if sensor == Sensors.nmeaLongitude:
                lon = read_nmea_coordinates(hex_segment[n : n + HEX_LEN_NMEA_LONGITUDE])
                results[HEX_TYPE_NMEA_LONGITUDE] = lon
                n += HEX_LEN_NMEA_LONGITUDE

            if sensor == Sensors.statusAndSign:
                signs = read_status_sign(hex_segment[n : n + HEX_LEN_NMEA_STATUS_AND_SIGN])
                if signs is np.nan:
                    results[HEX_TYPE_NMEA_LATITUDE] *= np.nan
                    results[HEX_TYPE_NMEA_LONGITUDE] *= np.nan
                else:
                    results[HEX_TYPE_NMEA_LATITUDE] *= signs[0]
                    results[HEX_TYPE_NMEA_LONGITUDE] *= signs[1]
                n += HEX_LEN_NMEA_STATUS_AND_SIGN

            if sensor == Sensors.nmeaTime:
                seconds_since_2000 = read_nmea_time(hex_segment[n : n + HEX_LEN_NMEA_TIME])
                # naive, no timezone conversion
                results[HEX_TYPE_NMEA_TIME] = datetime(1970, 1, 1) + timedelta(
                    seconds=seconds_since_2000 + SECONDS_BETWEEN_EPOCH_AND_2000
                )
                n += HEX_LEN_NMEA_TIME

    if moored_mode:
        seconds_since_2000 = int(hex_segment[n : n + HEX_LEN_DATE_TIME], 16)
        results[HEX_TYPE_DATE_TIME] = datetime.fromtimestamp(
            seconds_since_2000 + SECONDS_BETWEEN_EPOCH_AND_2000
        )
        n += HEX_LEN_DATE_TIME

    # Validate hex length. Ensure length matches what is expected based
    # on enabled sensors and moored mode. If it does not, set all values
    # to NaN
    
    if n != len(hex_segment.strip()):
        for key in results:
            results[key] = np.nan
        return results

    # Final validation to check if any values have been set to NaN, and if so,
    # set all values to NaN
    for key in results:
        if results[key] == np.nan:
            logger.warning("Invalid sample detected, values set to NaN")
            for key in results:
                results[key] = np.nan
            break
    return results


def read_SBE19plus_format_0(*args, **kwargs):
    warnings.warn("Deprecated, use read_sbe19plus_format_0", DeprecationWarning)
    return read_sbe19plus_format_0(*args, **kwargs)


def read_sbe37sm_format_0(
    hex_segment: str = "",
    enabled_sensors: Union[List[Sensors], None] = None,
) -> Dict[str, Union[int, float, datetime]]:
    """Converts a 37 family data hex string into engineering units.

    :param hex_segment: one line from a hex data file
    :param enabled_sensors: array of Sensors that are enabled. For 37
        this is always temperature, conductivity, pressure. Defaults to
        False

    :return: the 37 family sensor values in engineering units that were
        extracted from the input hex string
    """

    results: Dict[str, Union[int, float, datetime]] = {}
    n = 0
    results[HEX_TYPE_TEMPERATURE] = int(hex_segment[n:HEX_LEN_TEMPERATURE], 16)
    n += HEX_LEN_TEMPERATURE

    results[HEX_TYPE_CONDUCTIVITY] = int(hex_segment[n : n + HEX_LEN_CONDUCTIVITY], 16) / 256
    n += HEX_LEN_CONDUCTIVITY

    if enabled_sensors and Sensors.SBE63 in enabled_sensors:
        results[HEX_TYPE_SBE63_PHASE] = (
            int(hex_segment[n : n + HEX_LEN_SBE63_PHASE], 16) / 100000 - 10
        )
        n += HEX_LEN_SBE63_PHASE

        results[HEX_TYPE_SBE63_TEMPERATURE] = (
            int(hex_segment[n : n + HEX_LEN_TEMPERATURE], 16) / 1000000 - 1
        )
        n += HEX_LEN_TEMPERATURE

    if enabled_sensors and Sensors.Pressure in enabled_sensors:
        results[HEX_TYPE_PRESSURE] = int(hex_segment[n : n + HEX_LEN_PRESSURE], 16)
        n += HEX_LEN_PRESSURE

        result = int(hex_segment[n : n + HEX_LEN_TEMPERATURE_COMPENSATION], 16)
        results[HEX_TYPE_TEMPERATURE_COMPENSATION] = result
        n += HEX_LEN_TEMPERATURE_COMPENSATION

    seconds_since_2000 = int(hex_segment[n : n + HEX_LEN_DATE_TIME], 16)
    results[HEX_TYPE_DATE_TIME] = datetime(1970, 1, 1) + timedelta(
        seconds=seconds_since_2000 + SECONDS_BETWEEN_EPOCH_AND_2000
    )
    n += HEX_LEN_DATE_TIME

    # Validate hex length
    if n != len(hex_segment.strip()):
        raise ValueError("Hex string length does not match expectation based on enabled sensors")

    return results


def read_SBE37SM_format_0(*args, **kwargs):
    warnings.warn("Deprecated, use read_sbe37sm_format_0", DeprecationWarning)
    return read_sbe37sm_format_0(*args, **kwargs)


def read_nmea_coordinates(hex_segment: str):
    """Converts a 3 byte NMEA hex string to latitude or longitude

    :param hex_segment: 3 byte hex string
    :raises RuntimeWarning: raised if the hex string is the wrong length
    :return: latitude or longitide coordinate
    """
    if len(hex_segment) != 6:
        return np.nan

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
        return np.nan
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
        return [np.nan, np.nan]
    return signs


def read_nmea_time(hex_segment: str):
    """Convert an 8 byte hex string to the number of seconds since 2000

    :param hex_segment: an 8 byte hex string
    :raises RuntimeWarning: raised if the hex string is the wrong length
    :return: _description_
    """
    if len(hex_segment) != 8:
        return np.nan
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
    decimal_a = int(hex_segment[0:HEX_LEN_SBE911_VOLTAGE], 16)
    voltage_a = 5 * (1 - decimal_a / 4095)

    decimal_b = int(hex_segment[HEX_LEN_SBE911_VOLTAGE : HEX_LEN_SBE911_VOLTAGE * 2], 16)
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
