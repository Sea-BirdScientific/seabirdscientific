"""Functions for processing instrument data."""

import warnings
import xml.etree.ElementTree as ET
from datetime import datetime, timedelta
from enum import Enum
from itertools import dropwhile, takewhile
from logging import getLogger
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd
import xarray as xr

from seabirdscientific.constants import COUNTS_TO_VOLTS, SECONDS_BETWEEN_EPOCH_AND_2000
from seabirdscientific.utils import WarnAllMembersMeta

logger = getLogger(__name__)


class SensorID(Enum):
    ALTIMETER = 0
    CONDUCTIVITY = 3
    NOT_IN_USE = 27
    SBE43_DISSOLVED_OXYGEN = 38
    SBE27_ORP = 41
    SBE18_PH = 43
    DIGIQUARTZ_PRESSURE = 45
    STRAIN_GAUGE_PRESSURE = 46
    SPAR_BIOSPHERICAL = 51
    TEMPERATURE_FREQUENCY = 55
    SBE38_TEMPERATURE = 56
    TEMPERATURE = 58
    CUSTOM_POLYNOMIAL = 61
    OPTODE_OXYGEN = 64
    CSTAR = 71
    SBE63_DISSOLVED_OXYGEN_PHASE = 72
    PAR = 76
    SEAFET_PH_INTERNAL = 77
    SEAFET_PH_EXTERNAL = 78
    SPAR_LINEAR = 79
    CUSTOM_EXPONENTIAL = 80
    SPAR_LOGARITHMIC = 83
    SUNA = 86  # New sensors below here ----------------
    SEAFET_TEMPERATURE = 87
    INTERNAL_TEMPERATURE = 88
    RELATIVE_HUMIDITY = 89
    ECO_SERIAL_BACKSCATTERING_880 = 90
    ECO_SERIAL_BACKSCATTERING_870 = 91
    ECO_SERIAL_BACKSCATTERING_780 = 92
    ECO_SERIAL_BACKSCATTERING_770 = 93
    ECO_SERIAL_BACKSCATTERING_720 = 94
    ECO_SERIAL_BACKSCATTERING_715 = 95
    ECO_SERIAL_BACKSCATTERING_700 = 96
    ECO_SERIAL_BACKSCATTERING_676 = 97
    ECO_SERIAL_BACKSCATTERING_650 = 98
    ECO_SERIAL_BACKSCATTERING_630 = 99
    ECO_SERIAL_BACKSCATTERING_595 = 100
    ECO_SERIAL_BACKSCATTERING_532 = 101
    ECO_SERIAL_BACKSCATTERING_530 = 102
    ECO_SERIAL_BACKSCATTERING_510 = 103
    ECO_SERIAL_BACKSCATTERING_488 = 104
    ECO_SERIAL_BACKSCATTERING_470 = 105
    ECO_SERIAL_BACKSCATTERING_440 = 106
    ECO_SERIAL_BACKSCATTERING_412 = 107
    ECO_SERIAL_FDOM = 108
    ECO_SERIAL_FLUR = 109
    ECO_SERIAL_FLRH = 110
    ECO_SERIAL_FLPE = 111
    ECO_SERIAL_FLPC = 112
    ECO_SERIAL_CHL_435 = 113
    ECO_SERIAL_CHL_470 = 114
    CONDUCTIVITY_FREQUENCY = 115
    SBE27_PH = 116
    CUSTOM_EQUATION = 117
    PHYRE_PH_INTERNAL = 118
    PHYRE_PH_EXTERNAL = 119
    PHYRE_TEMPERATURE_INTERNAL = 120
    PHYRE_TEMPERATURE_EXTERNAL = 121
    PH_VK = 122
    PH_IB = 123
    PH_IK = 124
    STRAIN_GAUGE_PRESSURE_TEMPERATURE = 125
    DIGIQUARTZ_PRESSURE_TEMPERATURE = 126
    GTD_PRESSURE = 127  # TODO This may have an existing SensorID
    GTD_TEMPERATURE = 128  # TODO This may have an existing SensorID
    SBE63_DISSOLVED_OXYGEN_TEMPERATURE = 129
    NMEA_LATITUDE = 130
    NMEA_LONGITUDE = 131
    NMEA_STATUS_AND_SIGN = 132
    NMEA_DATE_TIME = 133
    INSTRUMENT_DATE_TIME = 134
    VOLTAGE_0 = 135  # raw voltage
    VOLTAGE_1 = 136
    VOLTAGE_2 = 137
    VOLTAGE_3 = 138
    VOLTAGE_4 = 139
    VOLTAGE_5 = 140
    VOLTAGE_6 = 141
    VOLTAGE_7 = 142
    ECO_VOLTAGE_BACKSCATTERING_880 = 143
    ECO_VOLTAGE_BACKSCATTERING_870 = 144
    ECO_VOLTAGE_BACKSCATTERING_780 = 145
    ECO_VOLTAGE_BACKSCATTERING_770 = 146
    ECO_VOLTAGE_BACKSCATTERING_720 = 147
    ECO_VOLTAGE_BACKSCATTERING_715 = 148
    ECO_VOLTAGE_BACKSCATTERING_700 = 149
    ECO_VOLTAGE_BACKSCATTERING_676 = 150
    ECO_VOLTAGE_BACKSCATTERING_650 = 151
    ECO_VOLTAGE_BACKSCATTERING_630 = 152
    ECO_VOLTAGE_BACKSCATTERING_595 = 153
    ECO_VOLTAGE_BACKSCATTERING_532 = 154
    ECO_VOLTAGE_BACKSCATTERING_530 = 155
    ECO_VOLTAGE_BACKSCATTERING_510 = 156
    ECO_VOLTAGE_BACKSCATTERING_488 = 157
    ECO_VOLTAGE_BACKSCATTERING_470 = 158
    ECO_VOLTAGE_BACKSCATTERING_440 = 159
    ECO_VOLTAGE_BACKSCATTERING_412 = 160
    ECO_VOLTAGE_FDOM = 161
    ECO_VOLTAGE_FLUR = 162
    ECO_VOLTAGE_FLRH = 163
    ECO_VOLTAGE_FLPE = 164
    ECO_VOLTAGE_FLPC = 165
    ECO_VOLTAGE_CHL_435 = 166
    ECO_VOLTAGE_CHL_470 = 167
    SBE911_PUMP_STATUS = 168
    SBE911_BOTTOM_CONTACT_STATUS = 169
    SBE911_CONFIRM_STATUS = 170
    SBE911_MODEM_STATUS = 171
    DATA_INTEGRITY = 172
    ERROR_FLAG = 173
    VOLTAGE = 998
    UNSUPORTED = 999


"""Possible lengths for hex data types"""
HEX_LEN_TEMPERATURE = 6
HEX_LEN_SBE38_TEMPERATURE = 6
HEX_LEN_CONDUCTIVITY = 6
HEX_LEN_STRAIN_GAUGE_PRESSURE = 6
HEX_LEN_STRAIN_GAUGE_PRESSURE_TEMP_COMP = 4
HEX_LEN_DIGIQUARTZ_PRESSURE = 6
HEX_LEN_DIGIQUARTZ_PRESSURE_TEMP_COMP = 3
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
HEX_LEN_SBE63_TEMPERATURE = 6
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
HEX_TYPE_PH_TEMPERATURE_EXTERNAL = "ph temperature external"
HEX_TYPE_PH_TEMPERATURE_INTERNAL = "ph temperature internal"
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
    pHyre = "pHyre"
    pHyreCAT = "pHyreCAT"
    HydroCAT = "HydroCAT"
    HydroCATODO = "HydroCAT-ODO"


class Sensors(Enum):
    """Available sensors to read hex data from"""

    Temperature = "Temperature"
    SecondaryTemperature = "SecondaryTemperature"
    Conductivity = "Conductivity"
    SecondaryConductivity = "SecondaryConductivity"
    Pressure = "Pressure"
    DigiquartzPressure = "DigiquartzPressure"
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


def _read_seasoft_cnv_file(filepath: Path | str) -> xr.Dataset:
    """Import the data from a .cnv file and put it into an xarray
    Dataset. Duplicate varioable names will have a number appended. For
    example, the second "depSM" becomes "depSM_1".

    :param filepath: the path to the .cnv file to be imported

    :return: the imported data from the .cnv file

    """

    dataset = xr.Dataset({}, attrs={"file_name": Path(filepath).name})

    total_scans = 0
    data_lines = []

    logger.info("Unpacking instrument data from file: %s", filepath)

    with open(filepath, mode="r") as cnv:
        for n_line, line in enumerate(cnv):
            if line.startswith("# nvalues = "):
                total_scans = int(line[line.find("= ") + 2 : line.find("\n")])

            elif line.startswith("# name "):
                name = line[line.find("= ") + 2 : line.find(":")]

                # scan is already added to each data array
                if name == "scan":
                    name = "scan_count"

                left_bracket = line.find(" [")
                if left_bracket > 0:
                    long_name = line[line.find(": ") + 2 : left_bracket]
                    units = line[line.find("[") + 1 : line.find("]")]
                else:
                    long_name = line[line.find(": ") + 2 : line.find("\n")].strip()
                    units = ""

                safe_name = name
                # check for label collision
                # if collision, increment until an available key is found
                n_name = 1
                while safe_name in list(dataset.data_vars):
                    safe_name = f"{name}_{n_name}"
                    n_name += 1

                if safe_name != name:
                    logger.warning(f'Duplicate measurand "{name}" will use key "{safe_name}"')

                data_array = xr.DataArray(
                    data=np.zeros(total_scans),
                    dims=["scan"],
                    coords={"scan": np.arange(total_scans)},
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
                dataset.attrs["start_time"] = start_time.isoformat()

            elif line.startswith("*END*"):
                data_lines = cnv.readlines()
                break

    np_data = np.array(
        [np.fromstring(dl.replace("-", " -").replace("e -", "e-"), sep=" ") for dl in data_lines]
    )

    for n, measurand in enumerate(list(dataset.data_vars)):
        dataset[measurand].data = np_data[:, n]

    return dataset


def _column_array(values: list | tuple) -> np.ndarray:
    """Converts one column of parsed values to a numpy array

    Numeric columns become float64. ValueError covers strings that do not
    parse as floats, TypeError covers values that are already Python objects
    (datetimes), and both fall back to an object array.

    :param values: every value for a single measurand, in scan order

    :return: a 1D array of the column
    """
    try:
        return np.array(values, dtype=float)
    except (ValueError, TypeError):
        return np.array(values, dtype=object)


def _read_fathom_cnv_file(filepath: Path | str) -> xr.Dataset:
    dataset = xr.Dataset({}, attrs={"file_name": Path(filepath).name})
    with open(filepath) as f:
        block = takewhile(
            lambda line: not line.startswith("*END*"),
            dropwhile(lambda line: not line.startswith("# <FathomProcessing"), f),
        )
        fathom_xml = ET.fromstring("".join(line.lstrip("#* ") for line in block))
        dc_element = fathom_xml.find("Modules").find("DataConversion")
        export_element = fathom_xml.find("Modules").find("Export")
        total_scans = int(dc_element.find("NumScans").text)
        dataset = dataset.assign_coords(scan=np.arange(total_scans))
        ds_attrs = {
            "sample_interval": float(dc_element.find("Interval").get("Value")),
            "start_time": dc_element.find("StartTime").get("DateTime"),
        }
        dataset.attrs.update(ds_attrs)

        rows = [line.split() for line in f if line.strip()]
        columns = list(export_element.find("Columns"))

        if len(rows) != total_scans or any(len(row) != len(columns) for row in rows):
            raise ValueError(
                f"expected {total_scans}x{len(columns)} data block, "
                f"got {len(rows)}x{len(rows[0]) if rows else 0}"
            )

        data = list(zip(*rows))

        for n, column in enumerate(columns):
            name = column.attrib["ID"]
            attrs = {
                "sbs_name": name,
                "long_name": column.find("Name").text,
                "units": column.find("Units").text or "",
            }
            values = data[n]

            # update the scan coord if there is a scan variable in the data
            if name in dataset.coords:
                dataset[name] = dataset[name].copy(data=np.array(values, dtype=np.int64))
                dataset[name].attrs.update(attrs)
                continue

            safe_name = name
            # check for label collision
            # if collision, increment until an available key is found
            n_name = 1
            while safe_name in dataset.data_vars:
                safe_name = f"{name}_{n_name}"
                n_name += 1

            if safe_name != name:
                logger.warning(f'Duplicate measurand "{name}" will use key "{safe_name}"')

            dataset[safe_name] = xr.DataArray(
                data=_column_array(values),
                dims=["scan"],
                attrs=attrs,
            )

    return dataset


def read_cnv_file(
    filepath: Path | str, software: Literal["fathom", "seasoft"] = "fathom"
) -> xr.Dataset:
    if software == "fathom":
        return _read_fathom_cnv_file(filepath)
    elif software == "seasoft":
        return _read_seasoft_cnv_file(filepath)
    else:
        raise ValueError(f"Unknown software: {software} (must be fathom or seasoft)")


def cnv_to_instrument_data(filepath: Path | str) -> xr.Dataset:
    warnings.warn("Deprecated, use read_cnv_file", DeprecationWarning)
    return read_cnv_file(filepath, "seasoft")


def read_hex_file(
    filepath: Path | str,
    instrument_type: InstrumentType,
    enabled_sensor_ids: list[SensorID],
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
                break

    # drop blank lines
    data_lines = [line for line in data_lines if not line.isspace()]

    attrs = {"file_name": Path(filepath).name}

    if not data_lines:
        return xr.Dataset({}, attrs=attrs)

    # Accumulate one list per measurand, then convert each column in a single
    # pass. Appending is cheaper than writing into a preallocated array cell by
    # cell, and avoids holding a scans x measurands array of boxed values.
    columns: dict[str, list] = {}

    for line in data_lines:
        hex_data = read_hex(
            line,
            instrument_type,
            enabled_sensor_ids,
        )
        if not columns:
            columns = {key: [] for key in hex_data}
        for key, column in columns.items():
            column.append(hex_data[key])

    # Build the Dataset in one construction so the scan coord is allocated once
    # and aligned once, rather than per measurand.
    return xr.Dataset(
        data_vars={
            key: xr.DataArray(
                data=_column_array(column),
                dims=["scan"],
                attrs={"sbs_name": key, "long_name": "", "units": ""},
            )
            for key, column in columns.items()
        },
        coords={"scan": np.arange(len(data_lines))},
        attrs=attrs,
    )


def read_hex(
    hex_segment: str, instrument_type: InstrumentType, enabled_sensor_ids: list[SensorID]
) -> dict[str, float | datetime]:
    """Converts a data hex string into raw decimal values"""

    results: dict[str, int | float | datetime] = {}
    n = 0
    sbe911_status: int | None = None
    replicate_by_sensor_id = [0] * 1000

    for sensor_id in enabled_sensor_ids:
        key = str(sensor_id) + "_" + str(replicate_by_sensor_id[int(sensor_id.value)])
        # -------- Internal Sensors  ---------
        if sensor_id in [
            SensorID.TEMPERATURE,
            SensorID.SEAFET_TEMPERATURE,
            SensorID.PHYRE_TEMPERATURE_EXTERNAL,
            SensorID.PHYRE_TEMPERATURE_INTERNAL,
        ]:
            if instrument_type == InstrumentType.SBE911Plus:
                results[key] = int(hex_segment[n : n + HEX_LEN_TEMPERATURE], 16) / 256
            else:
                results[key] = int(hex_segment[n : n + HEX_LEN_TEMPERATURE], 16)
            n += HEX_LEN_TEMPERATURE

        elif sensor_id == SensorID.CONDUCTIVITY:
            if instrument_type == InstrumentType.pHyreCAT:
                results[key] = int(hex_segment[n : n + HEX_LEN_CONDUCTIVITY], 16)
            else:
                results[key] = int(hex_segment[n : n + HEX_LEN_CONDUCTIVITY], 16) / 256
            n += HEX_LEN_CONDUCTIVITY

        elif sensor_id == SensorID.STRAIN_GAUGE_PRESSURE:
            results[key] = int(hex_segment[n : n + HEX_LEN_STRAIN_GAUGE_PRESSURE], 16)
            n += HEX_LEN_STRAIN_GAUGE_PRESSURE

        elif sensor_id == SensorID.STRAIN_GAUGE_PRESSURE_TEMPERATURE:
            pressure_temperature = int(
                hex_segment[n : n + HEX_LEN_STRAIN_GAUGE_PRESSURE_TEMP_COMP], 16
            )

            if instrument_type in [InstrumentType.SBE19Plus, InstrumentType.SBE16Plus]:
                pressure_temperature /= COUNTS_TO_VOLTS

            results[key] = pressure_temperature
            n += HEX_LEN_STRAIN_GAUGE_PRESSURE_TEMP_COMP

        elif sensor_id == SensorID.DIGIQUARTZ_PRESSURE:
            results[key] = int(hex_segment[n : n + HEX_LEN_DIGIQUARTZ_PRESSURE], 16) / 256
            n += HEX_LEN_DIGIQUARTZ_PRESSURE

        elif sensor_id == SensorID.DIGIQUARTZ_PRESSURE_TEMPERATURE:
            results[key] = int(hex_segment[n : n + HEX_LEN_DIGIQUARTZ_PRESSURE_TEMP_COMP], 16)
            if instrument_type == InstrumentType.SBE911Plus:
                n += HEX_LEN_SBE911_TEMPERATURE_COMPENSATION
            else:
                n += HEX_LEN_DIGIQUARTZ_PRESSURE_TEMP_COMP + 1  # 4th hex character is not used

        elif sensor_id in [SensorID.SEAFET_PH_EXTERNAL, SensorID.PHYRE_PH_EXTERNAL]:
            results[key] = int(hex_segment[n : n + HEX_LEN_VRS_EXTERNAL], 16)
            n += HEX_LEN_VRS_EXTERNAL

        elif sensor_id in [SensorID.SEAFET_PH_INTERNAL, SensorID.PHYRE_PH_INTERNAL]:
            results[key] = int(hex_segment[n : n + HEX_LEN_VRS_INTERNAL], 16)
            n += HEX_LEN_VRS_INTERNAL

        # -------------- Aux Serial Sensors -------------------

        elif sensor_id == SensorID.SBE38_TEMPERATURE:
            results[key] = int(hex_segment[n : n + HEX_LEN_SBE38_TEMPERATURE], 16) / 100000 - 10
            n += HEX_LEN_SBE38_TEMPERATURE

        elif sensor_id in [
            SensorID.ECO_SERIAL_BACKSCATTERING_880,
            SensorID.ECO_SERIAL_BACKSCATTERING_870,
            SensorID.ECO_SERIAL_BACKSCATTERING_780,
            SensorID.ECO_SERIAL_BACKSCATTERING_770,
            SensorID.ECO_SERIAL_BACKSCATTERING_720,
            SensorID.ECO_SERIAL_BACKSCATTERING_715,
            SensorID.ECO_SERIAL_BACKSCATTERING_700,
            SensorID.ECO_SERIAL_BACKSCATTERING_676,
            SensorID.ECO_SERIAL_BACKSCATTERING_650,
            SensorID.ECO_SERIAL_BACKSCATTERING_630,
            SensorID.ECO_SERIAL_BACKSCATTERING_595,
            SensorID.ECO_SERIAL_BACKSCATTERING_532,
            SensorID.ECO_SERIAL_BACKSCATTERING_530,
            SensorID.ECO_SERIAL_BACKSCATTERING_510,
            SensorID.ECO_SERIAL_BACKSCATTERING_488,
            SensorID.ECO_SERIAL_BACKSCATTERING_470,
            SensorID.ECO_SERIAL_BACKSCATTERING_440,
            SensorID.ECO_SERIAL_BACKSCATTERING_412,
            SensorID.ECO_SERIAL_FDOM,
            SensorID.ECO_SERIAL_FLUR,
            SensorID.ECO_SERIAL_FLRH,
            SensorID.ECO_SERIAL_FLPE,
            SensorID.ECO_SERIAL_FLPC,
            SensorID.ECO_SERIAL_CHL_435,
            SensorID.ECO_SERIAL_CHL_470,
        ]:
            results[key] = int(hex_segment[n : n + HEX_LEN_WETLABS_SINGLE_SENSOR], 16)
            n += HEX_LEN_WETLABS_SINGLE_SENSOR

        elif sensor_id == SensorID.GTD_PRESSURE:
            results[key] = int(hex_segment[n : n + HEX_LEN_GTD_PRESSURE], 16) / 10000
            n += HEX_LEN_GTD_PRESSURE

        elif sensor_id == SensorID.GTD_TEMPERATURE:
            results[key] = int(hex_segment[n : n + HEX_LEN_TEMPERATURE], 16) / 10000 - 10
            n += HEX_LEN_TEMPERATURE

        elif sensor_id == SensorID.OPTODE_OXYGEN:
            results[key] = int(hex_segment[n : n + HEX_LEN_OPTODE_OXYGEN], 16) / 10000 - 10
            n += HEX_LEN_OPTODE_OXYGEN

        elif sensor_id == SensorID.SBE63_DISSOLVED_OXYGEN_PHASE:
            if instrument_type == InstrumentType.pHyreCAT:
                results[key] = int(hex_segment[n : n + HEX_LEN_SBE63_PHASE], 16)
            else:
                results[key] = int(hex_segment[n : n + HEX_LEN_SBE63_PHASE], 16) / 100000 - 10
            n += HEX_LEN_SBE63_PHASE

        elif sensor_id == SensorID.SBE63_DISSOLVED_OXYGEN_TEMPERATURE:
            if instrument_type == InstrumentType.pHyreCAT:
                results[key] = int(hex_segment[n : n + HEX_LEN_SBE63_TEMPERATURE], 16)
            else:
                results[key] = (
                    int(hex_segment[n : n + HEX_LEN_SBE63_TEMPERATURE], 16) / 1000000 - 1
                )
            n += HEX_LEN_SBE63_TEMPERATURE

        # -------------- Aux Voltage Sensors -----------------
        elif sensor_id in [
            SensorID.VOLTAGE_0,
            SensorID.VOLTAGE_1,
            SensorID.VOLTAGE_2,
            SensorID.VOLTAGE_3,
            SensorID.VOLTAGE_4,
            SensorID.VOLTAGE_5,
            SensorID.VOLTAGE_6,
            SensorID.VOLTAGE_7,
            SensorID.SBE18_PH,
            SensorID.SBE27_PH,
            SensorID.SBE27_ORP,
            SensorID.SBE43_DISSOLVED_OXYGEN,
            SensorID.PAR,
            SensorID.SPAR_BIOSPHERICAL,
            SensorID.SPAR_LINEAR,
            SensorID.SPAR_LOGARITHMIC,
            SensorID.SUNA,
            SensorID.ALTIMETER,
            SensorID.CSTAR,
            SensorID.CUSTOM_POLYNOMIAL,
            SensorID.CUSTOM_EXPONENTIAL,
            SensorID.CUSTOM_EQUATION,
            SensorID.ECO_VOLTAGE_BACKSCATTERING_880,
            SensorID.ECO_VOLTAGE_BACKSCATTERING_870,
            SensorID.ECO_VOLTAGE_BACKSCATTERING_780,
            SensorID.ECO_VOLTAGE_BACKSCATTERING_770,
            SensorID.ECO_VOLTAGE_BACKSCATTERING_720,
            SensorID.ECO_VOLTAGE_BACKSCATTERING_715,
            SensorID.ECO_VOLTAGE_BACKSCATTERING_700,
            SensorID.ECO_VOLTAGE_BACKSCATTERING_676,
            SensorID.ECO_VOLTAGE_BACKSCATTERING_650,
            SensorID.ECO_VOLTAGE_BACKSCATTERING_630,
            SensorID.ECO_VOLTAGE_BACKSCATTERING_595,
            SensorID.ECO_VOLTAGE_BACKSCATTERING_532,
            SensorID.ECO_VOLTAGE_BACKSCATTERING_530,
            SensorID.ECO_VOLTAGE_BACKSCATTERING_510,
            SensorID.ECO_VOLTAGE_BACKSCATTERING_488,
            SensorID.ECO_VOLTAGE_BACKSCATTERING_470,
            SensorID.ECO_VOLTAGE_BACKSCATTERING_440,
            SensorID.ECO_VOLTAGE_BACKSCATTERING_412,
            SensorID.ECO_VOLTAGE_FDOM,
            SensorID.ECO_VOLTAGE_FLUR,
            SensorID.ECO_VOLTAGE_FLRH,
            SensorID.ECO_VOLTAGE_FLPE,
            SensorID.ECO_VOLTAGE_FLPC,
            SensorID.ECO_VOLTAGE_CHL_435,
            SensorID.ECO_VOLTAGE_CHL_470,
        ]:
            if instrument_type == InstrumentType.SBE911Plus and sensor_id in [
                SensorID.SPAR_BIOSPHERICAL,
                SensorID.SPAR_LINEAR,
                SensorID.SPAR_LOGARITHMIC,
            ]:
                n += 3  # unused bits
                results[key] = int(hex_segment[n : n + HEX_LEN_SBE911_SURFACE_PAR], 16) / 819
                n += HEX_LEN_SBE911_SURFACE_PAR
            elif instrument_type == InstrumentType.SBE911Plus:
                decimal_a = int(hex_segment[n : n + HEX_LEN_SBE911_VOLTAGE], 16)
                results[key] = 5 * (1 - decimal_a / 4095)
                n += HEX_LEN_SBE911_VOLTAGE
            else:
                results[key] = int(hex_segment[n : n + HEX_LEN_VOLTAGE], 16) / COUNTS_TO_VOLTS
                n += HEX_LEN_VOLTAGE

        # ------------- Additional Parameters ---------------
        elif sensor_id == SensorID.NMEA_LATITUDE:
            results[key] = read_nmea_coordinates(hex_segment[n : n + HEX_LEN_NMEA_LATITUDE])
            n += HEX_LEN_NMEA_LATITUDE

        elif sensor_id == SensorID.NMEA_LONGITUDE:
            results[key] = read_nmea_coordinates(hex_segment[n : n + HEX_LEN_NMEA_LONGITUDE])
            n += HEX_LEN_NMEA_LONGITUDE

        elif sensor_id == SensorID.NMEA_STATUS_AND_SIGN:
            signs = read_status_sign(hex_segment[n : n + HEX_LEN_NMEA_STATUS_AND_SIGN])
            if np.any(pd.isna(signs)):
                key = str(sensor_id) + "_" + str(replicate_by_sensor_id)
                results[
                    str(SensorID.NMEA_LATITUDE)
                    + "_"
                    + str(replicate_by_sensor_id[int(sensor_id.value)])
                ] = np.nan
                results[
                    str(SensorID.NMEA_LONGITUDE)
                    + "_"
                    + str(replicate_by_sensor_id[int(sensor_id.value)])
                ] = np.nan
            else:
                results[
                    str(SensorID.NMEA_LATITUDE)
                    + "_"
                    + str(replicate_by_sensor_id[int(sensor_id.value)])
                ] *= signs[0]
                results[
                    str(SensorID.NMEA_LONGITUDE)
                    + "_"
                    + str(replicate_by_sensor_id[int(sensor_id.value)])
                ] *= signs[1]
            n += HEX_LEN_NMEA_STATUS_AND_SIGN

        elif sensor_id == SensorID.NMEA_DATE_TIME:
            seconds_since_2000 = read_nmea_date_time(hex_segment[n : n + HEX_LEN_NMEA_TIME])
            # naive, no timezone conversion
            results[key] = datetime(1970, 1, 1) + timedelta(
                seconds=seconds_since_2000 + SECONDS_BETWEEN_EPOCH_AND_2000
            )
            n += HEX_LEN_NMEA_TIME

        elif sensor_id == SensorID.INSTRUMENT_DATE_TIME:
            date_time_hex = hex_segment[n : n + HEX_LEN_DATE_TIME]
            if instrument_type == InstrumentType.SBE911Plus:
                seconds_since_1970 = int(reverse_hex_bytes(date_time_hex), 16)
                results[key] = datetime(1970, 1, 1) + timedelta(seconds=seconds_since_1970)
            else:
                seconds_since_2000 = int(date_time_hex, 16)
                results[key] = datetime(1970, 1, 1) + timedelta(
                    seconds=seconds_since_2000 + SECONDS_BETWEEN_EPOCH_AND_2000
                )
            n += HEX_LEN_DATE_TIME

        elif sensor_id == SensorID.SBE911_PUMP_STATUS:
            if sbe911_status is None:
                sbe911_status = int(hex_segment[n : n + HEX_LEN_SBE911_STATUS], 16)
                n += HEX_LEN_SBE911_STATUS
            results[key] = sbe911_status & 1

        elif sensor_id == SensorID.SBE911_BOTTOM_CONTACT_STATUS:
            if sbe911_status is None:
                sbe911_status = int(hex_segment[n : n + HEX_LEN_SBE911_STATUS], 16)
                n += HEX_LEN_SBE911_STATUS
            results[key] = (sbe911_status >> 1) & 1

        elif sensor_id == SensorID.SBE911_CONFIRM_STATUS:
            if sbe911_status is None:
                sbe911_status = int(hex_segment[n : n + HEX_LEN_SBE911_STATUS], 16)
                n += HEX_LEN_SBE911_STATUS
            results[key] = (sbe911_status >> 2) & 1

        elif sensor_id == SensorID.SBE911_MODEM_STATUS:
            if sbe911_status is None:
                sbe911_status = int(hex_segment[n : n + HEX_LEN_SBE911_STATUS], 16)
                n += HEX_LEN_SBE911_STATUS
            results[key] = (sbe911_status >> 3) & 1

        elif sensor_id == SensorID.DATA_INTEGRITY:
            results[key] = int(hex_segment[n : n + HEX_LEN_SBE911_DATA_INTEGRITY], 16)
            n += HEX_LEN_SBE911_DATA_INTEGRITY

        elif sensor_id == SensorID.PH_VK:
            results[key] = int(hex_segment[n : n + HEX_LEN_VK], 16)
            n += HEX_LEN_VK

        elif sensor_id == SensorID.PH_IB:
            results[key] = int(hex_segment[n : n + HEX_LEN_IB], 16)
            n += HEX_LEN_IB

        elif sensor_id == SensorID.PH_IK:
            results[key] = int(hex_segment[n : n + HEX_LEN_IK], 16)
            n += HEX_LEN_IK

        elif sensor_id == SensorID.RELATIVE_HUMIDITY:
            results[key] = int(hex_segment[n : n + HEX_LEN_RELATIVE_HUMIDITY] + "0", 16)
            n += HEX_LEN_RELATIVE_HUMIDITY

        elif sensor_id == SensorID.INTERNAL_TEMPERATURE:
            results[key] = int(hex_segment[n : n + HEX_LEN_INTERNAL_TEMPERATURE] + "0", 16)
            n += HEX_LEN_INTERNAL_TEMPERATURE

        elif sensor_id == SensorID.ERROR_FLAG:
            results[key] = int(hex_segment[n : n + HEX_LEN_ERROR_FLAG], 16)
            n += HEX_LEN_ERROR_FLAG

        replicate_by_sensor_id[int(sensor_id.value)] += 1

    # Validate hex length. Ensure length matches what is expected based
    # on enabled sensors and moored mode. If it does not, set all values
    # to NaN
    if n != len(hex_segment.strip()):
        for key in results:
            results[key] = np.nan
        return results

    # Final validation to check if any values have been set to NaN, and if so,
    # set all values to NaN
    for key, value in results.items():
        if pd.isna(value):
            logger.warning("Invalid scan detected, values set to NaN")
            for key in results:
                results[key] = np.nan
            break

    return results


def read_nmea_coordinates(hex_segment: str):
    """Converts a 3 byte NMEA hex string to latitude or longitude

    :param hex_segment: 3 byte hex string
    :return: latitude or longitide coordinate or nan if the hex length
        is wrong
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
    :return: a list of two integers (1 or -1) or nan if the hex length
        is wrong
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


def read_nmea_date_time(hex_segment: str):
    """Convert an 8 byte hex string to the number of seconds since 2000

    :param hex_segment: an 8 byte hex string
    :return: an integer number of seconds or nan if the hex length was wrong
    """
    if len(hex_segment) != 8:
        return np.nan
    byte0 = hex_segment[0:2]
    byte1 = hex_segment[2:4]
    byte2 = hex_segment[4:6]
    byte3 = hex_segment[6:8]
    reformatted = int(byte3 + byte2 + byte1 + byte0, 16)
    return reformatted


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
