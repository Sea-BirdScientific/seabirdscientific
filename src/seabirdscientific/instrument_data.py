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
    ECO_SERIAL_UNKNOWN = 115
    CONDUCTIVITY_FREQUENCY = 116
    SBE27_PH = 117
    CUSTOM_EQUATION = 118
    PHYRE_PH_INTERNAL = 119
    PHYRE_PH_EXTERNAL = 120
    PHYRE_TEMPERATURE_INTERNAL = 121
    PHYRE_TEMPERATURE_EXTERNAL = 122
    PH_VK = 123
    PH_IB = 124
    PH_IK = 125
    STRAIN_GAUGE_PRESSURE_TEMPERATURE = 126
    DIGIQUARTZ_PRESSURE_TEMPERATURE = 127
    GTD_PRESSURE = 128  # TODO This may have an existing SensorID
    GTD_TEMPERATURE = 129  # TODO This may have an existing SensorID
    SBE63_DISSOLVED_OXYGEN_TEMPERATURE = 130
    NMEA_LATITUDE = 131
    NMEA_LONGITUDE = 132
    NMEA_STATUS_AND_SIGN = 133
    NMEA_DATE_TIME = 134
    INSTRUMENT_DATE_TIME = 135
    VOLTAGE_0 = 136  # raw voltage
    VOLTAGE_1 = 137
    VOLTAGE_2 = 138
    VOLTAGE_3 = 139
    VOLTAGE_4 = 140
    VOLTAGE_5 = 141
    VOLTAGE_6 = 142
    VOLTAGE_7 = 143
    ECO_VOLTAGE_BACKSCATTERING_880 = 144
    ECO_VOLTAGE_BACKSCATTERING_870 = 145
    ECO_VOLTAGE_BACKSCATTERING_780 = 146
    ECO_VOLTAGE_BACKSCATTERING_770 = 147
    ECO_VOLTAGE_BACKSCATTERING_720 = 148
    ECO_VOLTAGE_BACKSCATTERING_715 = 149
    ECO_VOLTAGE_BACKSCATTERING_700 = 150
    ECO_VOLTAGE_BACKSCATTERING_676 = 151
    ECO_VOLTAGE_BACKSCATTERING_650 = 152
    ECO_VOLTAGE_BACKSCATTERING_630 = 153
    ECO_VOLTAGE_BACKSCATTERING_595 = 154
    ECO_VOLTAGE_BACKSCATTERING_532 = 155
    ECO_VOLTAGE_BACKSCATTERING_530 = 156
    ECO_VOLTAGE_BACKSCATTERING_510 = 157
    ECO_VOLTAGE_BACKSCATTERING_488 = 158
    ECO_VOLTAGE_BACKSCATTERING_470 = 159
    ECO_VOLTAGE_BACKSCATTERING_440 = 160
    ECO_VOLTAGE_BACKSCATTERING_412 = 161
    ECO_VOLTAGE_FDOM = 162
    ECO_VOLTAGE_FLUR = 163
    ECO_VOLTAGE_FLRH = 164
    ECO_VOLTAGE_FLPE = 165
    ECO_VOLTAGE_FLPC = 166
    ECO_VOLTAGE_CHL_435 = 167
    ECO_VOLTAGE_CHL_470 = 168
    ECO_VOLTAGE_UNKNOWN = 169
    SBE911_DIAGNOSTICS = 170
    SBE911_PUMP_STATUS = 171
    SBE911_BOTTOM_CONTACT_STATUS = 172
    SBE911_CONFIRM_STATUS = 173
    SBE911_MODEM_STATUS = 174
    DATA_INTEGRITY = 175
    ERROR_FLAG = 176
    counts = 997
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


def sensor_id_to_hex_length(sensorID: SensorID, instrument_type: InstrumentType):
    general_map = {
        SensorID.TEMPERATURE: HEX_LEN_TEMPERATURE,
        SensorID.SEAFET_TEMPERATURE: HEX_LEN_TEMPERATURE,
        SensorID.PHYRE_TEMPERATURE_EXTERNAL: HEX_LEN_TEMPERATURE,
        SensorID.PHYRE_TEMPERATURE_INTERNAL: HEX_LEN_TEMPERATURE,
        SensorID.CONDUCTIVITY: HEX_LEN_CONDUCTIVITY,
        SensorID.STRAIN_GAUGE_PRESSURE: HEX_LEN_STRAIN_GAUGE_PRESSURE,
        SensorID.STRAIN_GAUGE_PRESSURE_TEMPERATURE: HEX_LEN_STRAIN_GAUGE_PRESSURE_TEMP_COMP,
        SensorID.DIGIQUARTZ_PRESSURE: HEX_LEN_DIGIQUARTZ_PRESSURE,
        SensorID.DIGIQUARTZ_PRESSURE_TEMPERATURE: HEX_LEN_DIGIQUARTZ_PRESSURE_TEMP_COMP,
        SensorID.SEAFET_PH_EXTERNAL: HEX_LEN_VRS_EXTERNAL,
        SensorID.PHYRE_PH_EXTERNAL: HEX_LEN_VRS_EXTERNAL,
        SensorID.SEAFET_PH_INTERNAL: HEX_LEN_VRS_INTERNAL,
        SensorID.PHYRE_PH_INTERNAL: HEX_LEN_VRS_INTERNAL,
        SensorID.SBE38_TEMPERATURE: HEX_LEN_SBE38_TEMPERATURE,
        SensorID.ECO_SERIAL_BACKSCATTERING_880: HEX_LEN_WETLABS_SINGLE_SENSOR,
        SensorID.ECO_SERIAL_BACKSCATTERING_870: HEX_LEN_WETLABS_SINGLE_SENSOR,
        SensorID.ECO_SERIAL_BACKSCATTERING_780: HEX_LEN_WETLABS_SINGLE_SENSOR,
        SensorID.ECO_SERIAL_BACKSCATTERING_770: HEX_LEN_WETLABS_SINGLE_SENSOR,
        SensorID.ECO_SERIAL_BACKSCATTERING_720: HEX_LEN_WETLABS_SINGLE_SENSOR,
        SensorID.ECO_SERIAL_BACKSCATTERING_715: HEX_LEN_WETLABS_SINGLE_SENSOR,
        SensorID.ECO_SERIAL_BACKSCATTERING_700: HEX_LEN_WETLABS_SINGLE_SENSOR,
        SensorID.ECO_SERIAL_BACKSCATTERING_676: HEX_LEN_WETLABS_SINGLE_SENSOR,
        SensorID.ECO_SERIAL_BACKSCATTERING_650: HEX_LEN_WETLABS_SINGLE_SENSOR,
        SensorID.ECO_SERIAL_BACKSCATTERING_630: HEX_LEN_WETLABS_SINGLE_SENSOR,
        SensorID.ECO_SERIAL_BACKSCATTERING_595: HEX_LEN_WETLABS_SINGLE_SENSOR,
        SensorID.ECO_SERIAL_BACKSCATTERING_532: HEX_LEN_WETLABS_SINGLE_SENSOR,
        SensorID.ECO_SERIAL_BACKSCATTERING_530: HEX_LEN_WETLABS_SINGLE_SENSOR,
        SensorID.ECO_SERIAL_BACKSCATTERING_510: HEX_LEN_WETLABS_SINGLE_SENSOR,
        SensorID.ECO_SERIAL_BACKSCATTERING_488: HEX_LEN_WETLABS_SINGLE_SENSOR,
        SensorID.ECO_SERIAL_BACKSCATTERING_470: HEX_LEN_WETLABS_SINGLE_SENSOR,
        SensorID.ECO_SERIAL_BACKSCATTERING_440: HEX_LEN_WETLABS_SINGLE_SENSOR,
        SensorID.ECO_SERIAL_BACKSCATTERING_412: HEX_LEN_WETLABS_SINGLE_SENSOR,
        SensorID.ECO_SERIAL_FDOM: HEX_LEN_WETLABS_SINGLE_SENSOR,
        SensorID.ECO_SERIAL_FLUR: HEX_LEN_WETLABS_SINGLE_SENSOR,
        SensorID.ECO_SERIAL_FLRH: HEX_LEN_WETLABS_SINGLE_SENSOR,
        SensorID.ECO_SERIAL_FLPE: HEX_LEN_WETLABS_SINGLE_SENSOR,
        SensorID.ECO_SERIAL_FLPC: HEX_LEN_WETLABS_SINGLE_SENSOR,
        SensorID.ECO_SERIAL_CHL_435: HEX_LEN_WETLABS_SINGLE_SENSOR,
        SensorID.ECO_SERIAL_CHL_470: HEX_LEN_WETLABS_SINGLE_SENSOR,
        SensorID.ECO_SERIAL_UNKNOWN: HEX_LEN_WETLABS_SINGLE_SENSOR,
        SensorID.GTD_PRESSURE: HEX_LEN_GTD_PRESSURE,
        SensorID.GTD_TEMPERATURE: HEX_LEN_TEMPERATURE,
        SensorID.OPTODE_OXYGEN: HEX_LEN_OPTODE_OXYGEN,
        SensorID.SBE63_DISSOLVED_OXYGEN_PHASE: HEX_LEN_SBE63_PHASE,
        SensorID.SBE63_DISSOLVED_OXYGEN_TEMPERATURE: HEX_LEN_SBE63_TEMPERATURE,
        SensorID.NMEA_LATITUDE: 2,  # 3 bytes, total of 6 hex chacters
        SensorID.NMEA_LONGITUDE: 2,  # 3 bytes, total of 6 hex chacters
        SensorID.NMEA_STATUS_AND_SIGN: HEX_LEN_NMEA_STATUS_AND_SIGN,
        SensorID.NMEA_DATE_TIME: HEX_LEN_NMEA_TIME,
        SensorID.INSTRUMENT_DATE_TIME: HEX_LEN_DATE_TIME,
        SensorID.SBE911_PUMP_STATUS: HEX_LEN_SBE911_STATUS,
        SensorID.SBE911_BOTTOM_CONTACT_STATUS: HEX_LEN_SBE911_STATUS,
        SensorID.SBE911_CONFIRM_STATUS: HEX_LEN_SBE911_STATUS,
        SensorID.SBE911_MODEM_STATUS: HEX_LEN_SBE911_STATUS,
        SensorID.SBE911_DIAGNOSTICS: HEX_LEN_SBE911_STATUS,
        SensorID.DATA_INTEGRITY: HEX_LEN_SBE911_DATA_INTEGRITY,
        SensorID.PH_VK: HEX_LEN_VK,
        SensorID.PH_IB: HEX_LEN_IB,
        SensorID.PH_IK: HEX_LEN_IK,
        SensorID.RELATIVE_HUMIDITY: HEX_LEN_RELATIVE_HUMIDITY,
        SensorID.INTERNAL_TEMPERATURE: HEX_LEN_INTERNAL_TEMPERATURE,
        SensorID.ERROR_FLAG: HEX_LEN_ERROR_FLAG,
    }

    if sensorID in [
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
        SensorID.ECO_VOLTAGE_UNKNOWN,
    ]:
        if instrument_type == InstrumentType.SBE911Plus:
            return HEX_LEN_SBE911_VOLTAGE
        else:
            return HEX_LEN_VOLTAGE
    else:
        return general_map[sensorID]


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
    values_array = np.asarray(values)
    if np.issubdtype(values_array.dtype, np.datetime64):
        return values_array

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

    # list of tuples matching SensorID to the key to be used in results dict like 'SensorID.TEMPERATURE_0'
    sensor_layout = _build_hex_sensor_layout(enabled_sensor_ids, instrument_type)

    columns = batch_convert_hex_to_dec(data_lines, sensor_layout)

    columns = _apply_transforms(columns, instrument_type, sensor_layout)

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


def _build_hex_sensor_layout(
    enabled_sensor_ids: list[SensorID], instrument_type: InstrumentType
) -> list[tuple[SensorID, str, int]]:
    """Precompute each sensor's output key, suffixed with its zero-based replicate number."""
    replicate_by_sensor_id: dict[SensorID, int] = {}
    sensor_layout: list[tuple[SensorID, str, int]] = []
    for sensor_id in enabled_sensor_ids:
        replicate = replicate_by_sensor_id.get(sensor_id, 0)
        hex_length = sensor_id_to_hex_length(sensor_id, instrument_type)

        if sensor_id in [SensorID.NMEA_LONGITUDE, SensorID.NMEA_LATITUDE]:
            sensor_layout.append((f"{sensor_id}_byte0", f"{sensor_id}_byte0", hex_length))
            sensor_layout.append((f"{sensor_id}_byte1", f"{sensor_id}_byte1", hex_length))
            sensor_layout.append((f"{sensor_id}_byte2", f"{sensor_id}_byte2", hex_length))
        elif sensor_id in [
            SensorID.SPAR_BIOSPHERICAL,
            SensorID.SPAR_LINEAR,
            SensorID.SPAR_LOGARITHMIC,
        ]:
            sensor_layout.append(
                (f"{sensor_id}_unused_bytes", f"{sensor_id}_unused_bytes", 3)
            )  # first three bytes not used in spar, this column is deleted later
            sensor_layout.append((sensor_id, f"{sensor_id}_{replicate}", hex_length))
        elif sensor_id in [
            SensorID.NMEA_STATUS_AND_SIGN,
            SensorID.DATA_INTEGRITY,
            SensorID.SBE911_DIAGNOSTICS,
            SensorID.INSTRUMENT_DATE_TIME,
        ]:
            # Don't add replicate suffix
            sensor_layout.append((sensor_id, str(sensor_id), hex_length))
        else:
            sensor_layout.append((sensor_id, f"{sensor_id}_{replicate}", hex_length))

        replicate_by_sensor_id[sensor_id] = replicate + 1
    return sensor_layout


def batch_convert_hex_to_dec(
    lines: list[str], sensor_layout: list[tuple[str, str, int]]
) -> dict[str, np.ndarray]:
    """
    Processes all lines at once using NumPy.
    Splits and converts fixed-length hex data to decimal columns at maximum speed.
    """
    # 1. Clean and convert your raw text lines into a NumPy string array
    # Stripping and casting to an 'S' (byte-string) type drops Python string overhead
    cleaned_lines = [line.strip() for line in lines]
    max_len = max(len(s) for s in cleaned_lines)
    arr = np.array(cleaned_lines, dtype=f"S{max_len}")

    # 2. View the data as raw byte characters (columns)
    # This turns a list of strings into a 2D matrix of individual characters
    char_matrix = arr.view(np.uint8).reshape(len(arr), -1)

    results = {}
    offset = 0

    # 3. Process each sensor column in bulk
    for _, key, hex_length in sensor_layout:
        # Slice out the character block for this specific sensor across ALL rows
        sensor_chars = char_matrix[:, offset : offset + hex_length]

        # Fast vector math to convert ASCII bytes to numeric values
        is_letter = sensor_chars > 57
        values = sensor_chars - 48
        values[is_letter] -= 7  # Adjusts 'A'-'F' or 'a'-'f' (assuming uniform casing)

        # Generate powers of 16 multipliers for this specific hex length
        multipliers = 16 ** np.arange(hex_length - 1, -1, -1)

        # Matrix multiplication computes the decimal equivalent instantly for all rows
        results[key] = np.dot(values, multipliers)

        offset += hex_length

    return results


def _apply_transforms(
    columns: dict[str, np.ndarray],
    instrument_type: InstrumentType,
    sensor_layout: list[tuple[SensorID, str, int]],
) -> dict[str, float | datetime]:
    """Converts a data hex string into raw decimal values.

    :param hex_segment: the hex string for one scan
    :param instrument_type: the instrument that generated the scan
    :param sensor_layout: ordered (sensor id, output key) pairs from
        _build_hex_sensor_layout, prepared once per sensor configuration
    :return: a dictionary of output key to raw value
    """

    for sensor_id, key, hex_length in sensor_layout:
        # -------- Internal Sensors  ---------
        if sensor_id == SensorID.TEMPERATURE and instrument_type == InstrumentType.SBE911Plus:
            columns[key] = columns[key] / 256

        elif sensor_id == SensorID.CONDUCTIVITY:
            if instrument_type != InstrumentType.pHyreCAT:
                columns[key] = columns[key] / 256

        elif sensor_id == SensorID.STRAIN_GAUGE_PRESSURE_TEMPERATURE:
            if instrument_type in [InstrumentType.SBE19Plus, InstrumentType.SBE16Plus]:
                columns[key] = columns[key] / COUNTS_TO_VOLTS

        elif sensor_id == SensorID.DIGIQUARTZ_PRESSURE:
            columns[key] = columns[key] / 256

        # -------------- Aux Serial Sensors -------------------

        elif sensor_id == SensorID.SBE38_TEMPERATURE:
            columns[key] = columns[key] / 100000 - 10

        elif sensor_id == SensorID.GTD_PRESSURE:
            columns[key] = columns[key] / 10000

        elif sensor_id in [SensorID.GTD_TEMPERATURE, SensorID.OPTODE_OXYGEN]:
            columns[key] = columns[key] / 10000 - 10

        elif sensor_id == SensorID.SBE63_DISSOLVED_OXYGEN_PHASE:
            if instrument_type != InstrumentType.pHyreCAT:
                columns[key] = columns[key] / 100000 - 10

        elif sensor_id == SensorID.SBE63_DISSOLVED_OXYGEN_TEMPERATURE:
            if instrument_type != InstrumentType.pHyreCAT:
                columns[key] = columns[key] / 1000000 - 1

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
                # n += 3  # unused bits
                columns[key] = columns[key] / 819
                columns.pop(f"{sensor_id}_unused_bytes")
            elif instrument_type == InstrumentType.SBE911Plus:
                decimal_a = columns[key]
                columns[key] = 5 * (1 - decimal_a / 4095)
            else:
                columns[key] = columns[key] / COUNTS_TO_VOLTS

        # ------------- Additional Parameters ---------------
        elif sensor_id == f"{SensorID.NMEA_LATITUDE}_byte0":
            columns[str(SensorID.NMEA_LATITUDE)] = read_nmea_coordinates(
                columns[f"{SensorID.NMEA_LATITUDE}_byte0"],
                columns[f"{SensorID.NMEA_LATITUDE}_byte1"],
                columns[f"{SensorID.NMEA_LATITUDE}_byte2"],
            )
            columns.pop(f"{SensorID.NMEA_LATITUDE}_byte0")
            columns.pop(f"{SensorID.NMEA_LATITUDE}_byte1")
            columns.pop(f"{SensorID.NMEA_LATITUDE}_byte2")

        elif sensor_id == f"{SensorID.NMEA_LONGITUDE}_byte0":
            columns[str(SensorID.NMEA_LONGITUDE)] = read_nmea_coordinates(
                columns[f"{SensorID.NMEA_LONGITUDE}_byte0"],
                columns[f"{SensorID.NMEA_LONGITUDE}_byte1"],
                columns[f"{SensorID.NMEA_LONGITUDE}_byte2"],
            )
            columns.pop(f"{SensorID.NMEA_LONGITUDE}_byte0")
            columns.pop(f"{SensorID.NMEA_LONGITUDE}_byte1")
            columns.pop(f"{SensorID.NMEA_LONGITUDE}_byte2")

        elif sensor_id == SensorID.NMEA_STATUS_AND_SIGN:
            signs = read_status_sign(columns[key])
            if np.any(pd.isna(signs)):
                columns[str(SensorID.NMEA_LATITUDE)] = np.nan
                columns[str(SensorID.NMEA_LONGITUDE)] = np.nan
            else:
                columns[str(SensorID.NMEA_LATITUDE)] *= signs[0]
                columns[str(SensorID.NMEA_LONGITUDE)] *= signs[1]
            columns.pop(key)

        elif sensor_id == SensorID.NMEA_DATE_TIME:
            seconds_since_2000 = columns[key]
            seconds_since_2000 = (
                ((seconds_since_2000 & 0x000000FF) << 24)
                | ((seconds_since_2000 & 0x0000FF00) << 8)
                | ((seconds_since_2000 & 0x00FF0000) >> 8)
                | ((seconds_since_2000 & 0xFF000000) >> 24)
            )
            columns[key] = np.datetime64("1970-01-01T00:00:00", "s") + (
                seconds_since_2000 + SECONDS_BETWEEN_EPOCH_AND_2000
            ).astype("timedelta64[s]")

        elif sensor_id == SensorID.INSTRUMENT_DATE_TIME:
            if instrument_type == InstrumentType.SBE911Plus:
                seconds_since_1970 = columns[key]
                seconds_since_1970 = (
                    ((seconds_since_1970 & 0x000000FF) << 24)
                    | ((seconds_since_1970 & 0x0000FF00) << 8)
                    | ((seconds_since_1970 & 0x00FF0000) >> 8)
                    | ((seconds_since_1970 & 0xFF000000) >> 24)
                )
                columns[key] = np.datetime64("1970-01-01T00:00:00", "s") + (
                    seconds_since_1970
                ).astype("timedelta64[s]")

            else:
                seconds_since_2000 = columns[key]
                columns[key] = np.datetime64("1970-01-01T00:00:00", "s") + (
                    seconds_since_2000 + SECONDS_BETWEEN_EPOCH_AND_2000
                ).astype("timedelta64[s]")

        elif sensor_id in [SensorID.RELATIVE_HUMIDITY, SensorID.INTERNAL_TEMPERATURE]:
            columns[key] = columns[key] << 4

        elif sensor_id == SensorID.SBE911_DIAGNOSTICS:
            columns[str(SensorID.SBE911_PUMP_STATUS)] = columns[key] & 1
            columns[str(SensorID.SBE911_BOTTOM_CONTACT_STATUS)] = (columns[key] >> 1) & 1
            columns[str(SensorID.SBE911_CONFIRM_STATUS)] = (columns[key] >> 2) & 1
            columns[str(SensorID.SBE911_MODEM_STATUS)] = (columns[key] >> 3) & 1
            columns.pop(key)

    # Validate hex length. Ensure length matches what is expected based
    # on enabled sensors and moored mode. If it does not, set all values
    # to NaN
    # if n != len(hex_segment.strip()):
    #     for key in results:
    #         results[key] = np.nan
    #     return results

    # Final validation to check if any values have been set to NaN, and if so,
    # set all values to NaN
    # for key, value in results.items():
    #     if pd.isna(value):
    #         logger.warning("Invalid scan detected, values set to NaN")
    #         for key in results:
    #             results[key] = np.nan
    #         break

    return columns


def read_nmea_coordinates(byte0: int, byte1: int, byte2: int):
    """Converts a 3 byte NMEA hex string to latitude or longitude

    :param hex_segment: 3 byte hex string
    :return: latitude or longitide coordinate or nan if the hex length
        is wrong
    """
    coordinate = (byte0 * 65536 + byte1 * 256 + byte2) / 50000
    return coordinate


def read_status_sign(integer: np.ndarray) -> list[np.ndarray]:
    """Converts a hex byte to the signs for NMEA latitude and longitude

    :param hex_segment: 1 byte hex string
    :return: latitude and longitude sign arrays containing 1 or -1
    """
    latitude_sign = np.where((integer & 0x80) == 0, 1, -1)
    longitude_sign = np.where((integer & 0x40) == 0, 1, -1)
    return [latitude_sign, longitude_sign]


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
