"""instrument data conversion unit tests."""

from datetime import datetime
from pathlib import Path
from typing import ClassVar

import pandas as pd
import pytest

import seabirdscientific.instrument_data as si

test_data = Path("./tests/resources/test-data")


class TestReadCnvFile:
    def test_read_cnv_file_seasoft(self):
        file_path = test_data / "SBE37SM.cnv"
        dataset = si.read_cnv_file(file_path, "seasoft")

        expected_variables = {
            "tv290C",
            "cond0S/m",
            "prdM",
            "prdE",
            "depSM",
            "depSF",
            "sal00",
            "flag",
        }
        assert expected_variables == set(dataset.data_vars)
        assert dataset.attrs["sample_interval"] == 120
        assert dataset["tv290C"].name == "tv290C"
        assert dataset["tv290C"].attrs["long_name"] == "Temperature"
        assert dataset["tv290C"].attrs["units"] == "ITS-90, deg C"
        assert dataset["tv290C"][0] == 20.8989
        assert dataset["tv290C"][-1] == -45.7713

    def test_duplicate_labels_seasoft(self):
        """
        Verify that duplicate labels are handled correctly.
        This .cnv file is from the Australian Institute of Marine Science
        """
        file_path = test_data / "WQR086CFACLWDB.cnv"
        dataset = si.read_cnv_file(file_path, "seasoft")

        expected_variables = {
            "prDM",
            "t090C",
            "c0S/m",
            "scan_count",
            "depSM",
            "timeS",
            "cpar",
            "spar",
            "turbWETntu0",
            "flECO-AFL",
            "par",
            "t190C",
            "c1S/m",
            "depSM_1",  # this duplicate "depSM" was a problem, increments to "depSM_1"
            "sal00",
            "nbin",
            "flag",
        }

        assert expected_variables == set(dataset.data_vars)

        # MeasurementSeries should have the original label
        assert dataset["depSM_1"].attrs["sbs_name"] == "depSM"

        first_values = [dataset[key].sel(scan=0).item() for key in list(dataset.data_vars)]
        first_data_line = "1.006 28.0831 5.811991 5248 0.998 218.616 6.9285e+01 3.9379e+02 11.8057 0.8223 2.7316e+02 28.0829 5.811250 1.000 36.2666 31 0.0000e+00"
        expected_first_values = [float(v) for v in first_data_line.split()]
        assert first_values == expected_first_values

    def test_read_cnv_fathom(self):
        filepath = test_data / "SBE19plus_01906398_2019_07_15_0033.cnv"
        dataset = si.read_cnv_file(filepath)
        expected_variables = {
            "c0S/m",
            "prdM",
            "sbeox0ML/L",
            "flag",
            "conservativeTemperature",
            "absoluteSalinity",
            "potentialDensity",
            "t090C",
            "sbeox0V",
            "pH",
            "timeS",
        }

        assert dataset.attrs.get("sample_interval") == 0.25
        assert set(dataset.data_vars) == expected_variables
        assert dataset["c0S/m"].attrs.get("sbs_name") == "c0S/m"
        assert dataset["c0S/m"].attrs.get("long_name") == "Conductivity"
        assert dataset["c0S/m"].attrs.get("units") == "S/m"
        assert dataset["c0S/m"][0] == 0.0896673177452
        assert dataset["c0S/m"][-1] == 1.1574933709080


class TestReadHex37smpOdo:
    @pytest.fixture
    def hex_data(self):
        return si.read_hex(
            # dataupload response
            "03AB8A0A5B492968C81899270877D308B5327396D4",
            si.InstrumentType.SBE37SMPODO,
            [
                si.SensorID.TEMPERATURE,
                si.SensorID.CONDUCTIVITY,
                si.SensorID.SBE63_DISSOLVED_OXYGEN_PHASE,
                si.SensorID.SBE63_DISSOLVED_OXYGEN_TEMPERATURE,
                si.SensorID.STRAIN_GAUGE_PRESSURE,
                si.SensorID.STRAIN_GAUGE_PRESSURE_TEMPERATURE,
                si.SensorID.INSTRUMENT_DATE_TIME,
            ],
        )

    print(hex_data)

    def test_read_37SMPODO_hex(self, hex_data):
        # Test case from instrument sample
        # Instrument tss response with output format 0 (raw decimal): 240522,  2651.285, 554963, 2229, 17.138, 0.612071,  27 Oct 2026, 17:07:32
        assert hex_data[str(si.SensorID.TEMPERATURE) + "_0"] == 240522
        assert round(hex_data[str(si.SensorID.CONDUCTIVITY) + "_0"], 3) == 2651.285
        assert (
            round(
                hex_data[str(si.SensorID.SBE63_DISSOLVED_OXYGEN_PHASE) + "_0"],
                5,
            )
            == 17.138
        )
        assert (
            round(hex_data[str(si.SensorID.SBE63_DISSOLVED_OXYGEN_TEMPERATURE) + "_0"], 6)
            == 0.612071
        )
        assert hex_data[str(si.SensorID.STRAIN_GAUGE_PRESSURE) + "_0"] == 554963
        assert hex_data[str(si.SensorID.STRAIN_GAUGE_PRESSURE_TEMPERATURE) + "_0"] == 2229


class TestReadHex19plus:
    @pytest.fixture
    def hex_data(self):
        return si.read_hex(
            #  T      C      P      PT   V0   V1   V2   V3   V4   V5   S0   S1   S2   DateTime
            # "03EC41 0A918B 0826F9 68C8 08C2 4891 0001 0000 0003 0000 00BF 0032 1022"          Fathom upload - missing date time?
            # "03EC41 0A918B 0826F9 68C8 08C1 4890 0000 0000 0002 0000 00BF 0032 1022 32739C93" outputformat=0 (raw hex) getsamples:0,0
            "03EC410A918B0826F968C808C24891000100000003000000BF00321022",
            si.InstrumentType.SBE19Plus,
            [
                si.SensorID.TEMPERATURE,
                si.SensorID.CONDUCTIVITY,
                si.SensorID.STRAIN_GAUGE_PRESSURE,
                si.SensorID.STRAIN_GAUGE_PRESSURE_TEMPERATURE,
                si.SensorID.VOLTAGE_0,
                si.SensorID.VOLTAGE_1,
                si.SensorID.VOLTAGE_2,
                si.SensorID.VOLTAGE_3,
                si.SensorID.VOLTAGE_4,
                si.SensorID.VOLTAGE_5,
                si.SensorID.ECO_SERIAL_BACKSCATTERING_412,
                si.SensorID.ECO_SERIAL_BACKSCATTERING_440,
                si.SensorID.ECO_SERIAL_BACKSCATTERING_470,
                # si.SensorID.INSTRUMENT_DATE_TIME,
            ],
        )

    def test_read_19plus_hex(self, hex_data):
        # outputformat=2 (raw decimal) getsamples:0,0 response:  257089, 2705.543, 534265, 2.0465, 0.1711, 1.4173, 0.0001, 0.0000, 0.0002, 0.0000, 191, 50, 4130, 27 Oct 2026 17:32:03
        assert hex_data[str(si.SensorID.TEMPERATURE) + "_0"] == 257089
        assert round(hex_data[str(si.SensorID.CONDUCTIVITY) + "_0"], 3) == 2705.543
        assert hex_data[str(si.SensorID.STRAIN_GAUGE_PRESSURE) + "_0"] == 534265
        assert (
            round(hex_data[str(si.SensorID.STRAIN_GAUGE_PRESSURE_TEMPERATURE) + "_0"], 4) == 2.0465
        )
        assert round(hex_data[str(si.SensorID.VOLTAGE_0) + "_0"], 4) == 0.1711
        assert round(hex_data[str(si.SensorID.VOLTAGE_1) + "_0"], 4) == 1.4173
        assert round(hex_data[str(si.SensorID.VOLTAGE_2) + "_0"], 4) == 0.0001
        assert round(hex_data[str(si.SensorID.VOLTAGE_3) + "_0"], 4) == 0.0000
        assert round(hex_data[str(si.SensorID.VOLTAGE_4) + "_0"], 4) == 0.0002
        assert round(hex_data[str(si.SensorID.VOLTAGE_5) + "_0"], 4) == 0.0000
        assert hex_data[str(si.SensorID.ECO_SERIAL_BACKSCATTERING_412) + "_0"] == 191
        assert hex_data[str(si.SensorID.ECO_SERIAL_BACKSCATTERING_440) + "_0"] == 50
        assert hex_data[str(si.SensorID.ECO_SERIAL_BACKSCATTERING_470) + "_0"] == 4130
        # assert pd.to_datetime(hex_data[str(si.SensorID.INSTRUMENT_DATE_TIME) + "_0"]) == datetime(
        #     2026, 10, 27, 17, 32, 3
        # )

    def test_read_19plus_nmea_hex(self):

        hex_data = si.read_hex(
            "08E1A913C11D08085550CB9D40033C2518A45E391940D767A82B",  # first scan in 19plus_V2_CTD_nmea_data.hex
            si.InstrumentType.SBE19Plus,
            [
                si.SensorID.TEMPERATURE,
                si.SensorID.CONDUCTIVITY,
                si.SensorID.STRAIN_GAUGE_PRESSURE,
                si.SensorID.STRAIN_GAUGE_PRESSURE_TEMPERATURE,
                si.SensorID.VOLTAGE_0,
                si.SensorID.VOLTAGE_1,
                si.SensorID.NMEA_LATITUDE,
                si.SensorID.NMEA_LONGITUDE,
                si.SensorID.NMEA_STATUS_AND_SIGN,
                si.SensorID.NMEA_DATE_TIME,
            ],
        )
        assert hex_data[str(si.SensorID.NMEA_LATITUDE) + "_0"] == 48.62280000
        assert hex_data[str(si.SensorID.NMEA_LONGITUDE) + "_0"] == -123.50002000
        assert pd.to_datetime(hex_data[str(si.SensorID.NMEA_DATE_TIME) + "_0"]) == datetime(
            2023, 3, 18, 12, 15, 19
        )


class TestReadHex16plusDigiquartz:
    def test_read_16plus_hex_digiquartz(self):

        hex_data = si.read_hex(
            "051E3C0A69C708726649CD0003000300040004000900040000640003E83273CF39",  # Fathom upload
            si.InstrumentType.SBE16Plus,
            [
                si.SensorID.TEMPERATURE,
                si.SensorID.CONDUCTIVITY,
                si.SensorID.STRAIN_GAUGE_PRESSURE,
                si.SensorID.STRAIN_GAUGE_PRESSURE_TEMPERATURE,
                si.SensorID.VOLTAGE_0,
                si.SensorID.VOLTAGE_1,
                si.SensorID.VOLTAGE_2,
                si.SensorID.VOLTAGE_3,
                si.SensorID.VOLTAGE_4,
                si.SensorID.VOLTAGE_5,
                si.SensorID.SBE63_DISSOLVED_OXYGEN_PHASE,
                si.SensorID.SBE63_DISSOLVED_OXYGEN_TEMPERATURE,
                si.SensorID.INSTRUMENT_DATE_TIME,
            ],
        )
        # outputformat=2 (raw decimal) getsamples:0,0 response: 335420, 2665.777, 553574, 1.4414, 0.0002, 0.0002, 0.0003, 0.0003, 0.0007, 0.0003, -9.999, -0.999000, 27 Oct 2026 21:08:09
        assert hex_data[str(si.SensorID.TEMPERATURE) + "_0"] == 335420
        assert round(hex_data[str(si.SensorID.CONDUCTIVITY) + "_0"], 3) == 2665.777
        assert hex_data[str(si.SensorID.STRAIN_GAUGE_PRESSURE) + "_0"] == 553574
        assert (
            round(hex_data[str(si.SensorID.STRAIN_GAUGE_PRESSURE_TEMPERATURE) + "_0"], 4) == 1.4414
        )
        assert round(hex_data[str(si.SensorID.VOLTAGE_0) + "_0"], 4) == 0.0002
        assert round(hex_data[str(si.SensorID.VOLTAGE_1) + "_0"], 4) == 0.0002
        assert round(hex_data[str(si.SensorID.VOLTAGE_2) + "_0"], 4) == 0.0003
        assert round(hex_data[str(si.SensorID.VOLTAGE_3) + "_0"], 4) == 0.0003
        assert round(hex_data[str(si.SensorID.VOLTAGE_4) + "_0"], 4) == 0.0007
        assert round(hex_data[str(si.SensorID.VOLTAGE_5) + "_0"], 4) == 0.0003
        assert hex_data[str(si.SensorID.SBE63_DISSOLVED_OXYGEN_PHASE) + "_0"] == -9.999
        assert hex_data[str(si.SensorID.SBE63_DISSOLVED_OXYGEN_TEMPERATURE) + "_0"] == -0.999000
        assert pd.to_datetime(hex_data[str(si.SensorID.INSTRUMENT_DATE_TIME) + "_0"]) == datetime(
            2026, 10, 27, 21, 8, 9
        )


class TestReadHex39plus:
    def test_read_39plus_hex(self):
        hex_data = si.read_hex(
            "32407EF51A485301A5418459",  # Fathom upload
            si.InstrumentType.SBE39Plus,
            [
                si.SensorID.INSTRUMENT_DATE_TIME,
                si.SensorID.TEMPERATURE,
                si.SensorID.STRAIN_GAUGE_PRESSURE,
                si.SensorID.STRAIN_GAUGE_PRESSURE_TEMPERATURE,
            ],
        )

        # outputformat=0 (raw decimal) getsamples:0,0 response: 1722451,  107841, 33881, 18 Sep 2026, 23:00:05
        assert hex_data[str(si.SensorID.TEMPERATURE) + "_0"] == 1722451
        assert hex_data[str(si.SensorID.STRAIN_GAUGE_PRESSURE) + "_0"] == 107841
        assert hex_data[str(si.SensorID.STRAIN_GAUGE_PRESSURE_TEMPERATURE) + "_0"] == 33881
        assert pd.to_datetime(hex_data[str(si.SensorID.INSTRUMENT_DATE_TIME) + "_0"]) == datetime(
            2026, 9, 18, 23, 0, 5
        )


class TestReadHexSBE911plus:
    # Expected values are from SBE Data Processing
    def test_read_911plus_hex(self):
        hex_data = si.read_hex(
            "122B2B19EEFC821642122BC818567B40DF5203788580976EF847110004941A26795B9280409C2353390B8868",
            si.InstrumentType.SBE911Plus,
            [
                si.SensorID.TEMPERATURE,
                si.SensorID.CONDUCTIVITY,
                si.SensorID.DIGIQUARTZ_PRESSURE,
                si.SensorID.TEMPERATURE,
                si.SensorID.CONDUCTIVITY,
                si.SensorID.VOLTAGE_0,
                si.SensorID.VOLTAGE_1,
                si.SensorID.VOLTAGE_2,
                si.SensorID.VOLTAGE_3,
                si.SensorID.VOLTAGE_4,
                si.SensorID.VOLTAGE_5,
                si.SensorID.VOLTAGE_6,
                si.SensorID.VOLTAGE_7,
                si.SensorID.SPAR_LINEAR,
                si.SensorID.NMEA_LATITUDE,
                si.SensorID.NMEA_LONGITUDE,
                si.SensorID.NMEA_STATUS_AND_SIGN,
                si.SensorID.DIGIQUARTZ_PRESSURE_TEMPERATURE,
                si.SensorID.SBE911_PUMP_STATUS,
                si.SensorID.SBE911_BOTTOM_CONTACT_STATUS,
                si.SensorID.SBE911_CONFIRM_STATUS,
                si.SensorID.SBE911_MODEM_STATUS,
                si.SensorID.DATA_INTEGRITY,
                si.SensorID.INSTRUMENT_DATE_TIME,
            ],
        )

        print(hex_data)

        assert round(hex_data[str(si.SensorID.TEMPERATURE) + "_0"], 3) == 4651.168
        assert round(hex_data[str(si.SensorID.CONDUCTIVITY) + "_0"], 3) == 6638.984
        assert round(hex_data[str(si.SensorID.DIGIQUARTZ_PRESSURE) + "_0"], 3) == 33302.258
        assert round(hex_data[str(si.SensorID.TEMPERATURE) + "_1"], 3) == 4651.781
        assert round(hex_data[str(si.SensorID.CONDUCTIVITY) + "_1"], 3) == 6230.480
        assert round(hex_data[str(si.SensorID.VOLTAGE_0) + "_0"], 4) == 3.7338
        assert round(hex_data[str(si.SensorID.VOLTAGE_1) + "_0"], 4) == 0.2112
        assert round(hex_data[str(si.SensorID.VOLTAGE_2) + "_0"], 4) == 4.9328
        assert round(hex_data[str(si.SensorID.VOLTAGE_3) + "_0"], 4) == 2.3370
        assert round(hex_data[str(si.SensorID.VOLTAGE_4) + "_0"], 4) == 2.4884
        assert round(hex_data[str(si.SensorID.VOLTAGE_5) + "_0"], 4) == 2.6777
        assert round(hex_data[str(si.SensorID.VOLTAGE_6) + "_0"], 4) == 0.1502
        assert round(hex_data[str(si.SensorID.VOLTAGE_7) + "_0"], 4) == 2.7912
        assert round(hex_data[str(si.SensorID.SPAR_LINEAR) + "_0"], 5) == 1.43101
        assert hex_data[str(si.SensorID.NMEA_LATITUDE) + "_0"] == 34.2757
        assert hex_data[str(si.SensorID.NMEA_LONGITUDE) + "_0"] == -120.0256
        # Status nibble 0b0011: pump on, no bottom contact (pump and bottom contact match SBE
        # Data Processing; Fathom read the bits most significant first)
        assert int(hex_data[str(si.SensorID.SBE911_PUMP_STATUS) + "_0"]) == 1
        assert int(hex_data[str(si.SensorID.SBE911_BOTTOM_CONTACT_STATUS) + "_0"]) == 1
        assert int(hex_data[str(si.SensorID.SBE911_CONFIRM_STATUS) + "_0"]) == 0
        assert int(hex_data[str(si.SensorID.SBE911_MODEM_STATUS) + "_0"]) == 0
        assert int(hex_data[str(si.SensorID.DATA_INTEGRITY) + "_0"]) == 83
        assert pd.to_datetime(hex_data[str(si.SensorID.INSTRUMENT_DATE_TIME) + "_0"]) == datetime(
            2025, 7, 28, 23, 43, 53
        )


class TestReadHexSBE911PlusStatusBits:
    # 5 frequencies, 4 voltage words, pressure sensor temperature, status, modulo count
    enabled_sensor_ids: ClassVar[list[si.SensorID]] = [
        si.SensorID.TEMPERATURE,
        si.SensorID.CONDUCTIVITY,
        si.SensorID.DIGIQUARTZ_PRESSURE,
        si.SensorID.TEMPERATURE,
        si.SensorID.CONDUCTIVITY,
        si.SensorID.VOLTAGE_0,
        si.SensorID.VOLTAGE_1,
        si.SensorID.VOLTAGE_2,
        si.SensorID.VOLTAGE_3,
        si.SensorID.VOLTAGE_4,
        si.SensorID.VOLTAGE_5,
        si.SensorID.VOLTAGE_6,
        si.SensorID.VOLTAGE_7,
        si.SensorID.DIGIQUARTZ_PRESSURE_TEMPERATURE,
        si.SensorID.SBE911_PUMP_STATUS,
        si.SensorID.SBE911_BOTTOM_CONTACT_STATUS,
        si.SensorID.SBE911_CONFIRM_STATUS,
        si.SensorID.SBE911_MODEM_STATUS,
        si.SensorID.DATA_INTEGRITY,
    ]

    def scan(self, nibble: int) -> str:
        return "0" * 30 + "0" * 24 + "ABC" + f"{nibble:X}" + "6C"

    @pytest.mark.parametrize("nibble", range(16))
    def test_bits_read_least_significant_first(self, nibble):
        result = si.read_hex(
            self.scan(nibble), si.InstrumentType.SBE911Plus, self.enabled_sensor_ids
        )
        assert result[str(si.SensorID.SBE911_PUMP_STATUS) + "_0"] == nibble & 1
        assert result[str(si.SensorID.SBE911_BOTTOM_CONTACT_STATUS) + "_0"] == (nibble >> 1) & 1
        assert result[str(si.SensorID.SBE911_CONFIRM_STATUS) + "_0"] == (nibble >> 2) & 1
        assert result[str(si.SensorID.SBE911_MODEM_STATUS) + "_0"] == (nibble >> 3) & 1
        assert result[str(si.SensorID.DIGIQUARTZ_PRESSURE_TEMPERATURE) + "_0"] == 0xABC
        assert result[str(si.SensorID.DATA_INTEGRITY) + "_0"] == 0x6C

    @pytest.mark.parametrize(
        "nibble, pump, confirm",
        [
            (0b0010, 0, 0),  # deck test: pump off
            (0b0011, 1, 0),  # in the water: pump on
            (0b0111, 1, 1),  # bottle firing
        ],
    )
    def test_status_seen_in_real_casts(self, nibble, pump, confirm):
        result = si.read_hex(
            self.scan(nibble), si.InstrumentType.SBE911Plus, self.enabled_sensor_ids
        )
        assert result[str(si.SensorID.SBE911_PUMP_STATUS) + "_0"] == pump
        assert result[str(si.SensorID.SBE911_CONFIRM_STATUS) + "_0"] == confirm
        # no bottom contact, modem carrier detected
        assert result[str(si.SensorID.SBE911_BOTTOM_CONTACT_STATUS) + "_0"] == 1
        assert result[str(si.SensorID.SBE911_MODEM_STATUS) + "_0"] == 0


class TestReadHexSeaFET2:
    def test_read_shallow_seaFET2_hex(self):
        hex_data = si.read_hex(
            "470F0245403C03884649C5CC7FD8357FF109FFF5D8327403BA0000",
            si.InstrumentType.SeaFET2,
            [
                si.SensorID.SEAFET_PH_EXTERNAL,
                si.SensorID.SEAFET_PH_INTERNAL,
                si.SensorID.SEAFET_TEMPERATURE,
                si.SensorID.PH_VK,
                si.SensorID.PH_IB,
                si.SensorID.PH_IK,
                si.SensorID.RELATIVE_HUMIDITY,
                si.SensorID.INTERNAL_TEMPERATURE,
                si.SensorID.INSTRUMENT_DATE_TIME,
                si.SensorID.ERROR_FLAG,
            ],
        )

        # outputformat=0 (raw decimal) getsamples:1,1 response: SEAFET00391,2026-10-28T00:52:10, 0000, 4656898, 4538428,  231494, 4834764, 8378421, 8384777, 23936, 65520
        assert hex_data[str(si.SensorID.SEAFET_PH_EXTERNAL) + "_0"] == 4656898
        assert hex_data[str(si.SensorID.SEAFET_PH_INTERNAL) + "_0"] == 4538428
        assert hex_data[str(si.SensorID.SEAFET_TEMPERATURE) + "_0"] == 231494
        assert hex_data[str(si.SensorID.PH_VK) + "_0"] == 4834764
        assert hex_data[str(si.SensorID.PH_IB) + "_0"] == 8378421
        assert hex_data[str(si.SensorID.PH_IK) + "_0"] == 8384777
        assert hex_data[str(si.SensorID.RELATIVE_HUMIDITY) + "_0"] == 65520
        assert hex_data[str(si.SensorID.INTERNAL_TEMPERATURE) + "_0"] == 23936
        assert hex_data[str(si.SensorID.INSTRUMENT_DATE_TIME) + "_0"] == datetime(
            2026, 10, 28, 0, 52, 10
        )
        assert hex_data[str(si.SensorID.ERROR_FLAG) + "_0"] == 0000


class TestReadHexSeapHOx:
    def test_read_shallow_seaphox_hex(self):
        hex_data = si.read_hex(
            "0396ED0A3EB708738908922963501869FC470F8D4541F4038B7049C8187FD7B17FF088FFF5DA32740BB30000",
            si.InstrumentType.SeapHox2,
            [
                si.SensorID.TEMPERATURE,
                si.SensorID.CONDUCTIVITY,
                si.SensorID.STRAIN_GAUGE_PRESSURE,
                si.SensorID.STRAIN_GAUGE_PRESSURE_TEMPERATURE,
                si.SensorID.SBE63_DISSOLVED_OXYGEN_PHASE,
                si.SensorID.SBE63_DISSOLVED_OXYGEN_TEMPERATURE,
                si.SensorID.SEAFET_PH_EXTERNAL,
                si.SensorID.SEAFET_PH_INTERNAL,
                si.SensorID.SEAFET_TEMPERATURE,
                si.SensorID.PH_VK,
                si.SensorID.PH_IB,
                si.SensorID.PH_IK,
                si.SensorID.RELATIVE_HUMIDITY,
                si.SensorID.INTERNAL_TEMPERATURE,
                si.SensorID.INSTRUMENT_DATE_TIME,
                si.SensorID.ERROR_FLAG,
            ],
        )
        #                                                       frame       date time            error  T       vrs_ext  vrs_int   ph_temp vk       ib       ik        p        pt    c         ox_ph   ox_t      it     rh
        # outputformat=0 (raw decimal) getsamples:1,1 response: SSPHOX00391,2026-10-28T01:26:11, 0000,  235245, 4657037, 4538868,  232304, 4835352, 8378289, 8384648,  553865,  2194, 2622.715, 17.124, 0.599996, 23968, 65520
        assert hex_data[str(si.SensorID.TEMPERATURE) + "_0"] == 235245
        assert round(hex_data[str(si.SensorID.CONDUCTIVITY) + "_0"], 3) == 2622.715
        assert hex_data[str(si.SensorID.STRAIN_GAUGE_PRESSURE) + "_0"] == 553865
        assert hex_data[str(si.SensorID.STRAIN_GAUGE_PRESSURE_TEMPERATURE) + "_0"] == 2194
        assert round(hex_data[str(si.SensorID.SBE63_DISSOLVED_OXYGEN_PHASE) + "_0"], 3) == 17.124
        assert (
            round(hex_data[str(si.SensorID.SBE63_DISSOLVED_OXYGEN_TEMPERATURE) + "_0"], 6)
            == 0.599996
        )

        assert hex_data[str(si.SensorID.SEAFET_PH_EXTERNAL) + "_0"] == 4657037
        assert hex_data[str(si.SensorID.SEAFET_PH_INTERNAL) + "_0"] == 4538868
        assert hex_data[str(si.SensorID.SEAFET_TEMPERATURE) + "_0"] == 232304
        assert hex_data[str(si.SensorID.PH_VK) + "_0"] == 4835352
        assert hex_data[str(si.SensorID.PH_IB) + "_0"] == 8378289
        assert hex_data[str(si.SensorID.PH_IK) + "_0"] == 8384648
        assert hex_data[str(si.SensorID.RELATIVE_HUMIDITY) + "_0"] == 65520
        assert hex_data[str(si.SensorID.INTERNAL_TEMPERATURE) + "_0"] == 23968
        assert hex_data[str(si.SensorID.INSTRUMENT_DATE_TIME) + "_0"] == datetime(
            2026, 10, 28, 1, 26, 11
        )
        assert hex_data[str(si.SensorID.ERROR_FLAG) + "_0"] == 0000


class TestReadpHyre:
    # The following expected values were output directly by a pHyre with the TSS command and
    # outputformat=0. The samples were then uploaded using the uploaddata command in terminal.
    @pytest.mark.parametrize(
        "hex_segment, date, event, ph_temperature_external, ph_temperature_internal, "
        "vrs_external, vrs_internal, vk, ib, ik",
        [
            (
                "4FB8A04FB3290AA3EE0AA42258D2A357994E99004B326176670008",
                "2026-10-13T23:08:23",
                8,
                697326,
                697378,
                5224608,
                5223209,
                5821091,
                5740878,
                10027083,
            ),
            (
                "4FB8A24FB3250AA3D50AA42258D28D57986C98FD553261766B0008",
                "2026-10-13T23:08:27",
                8,
                697301,
                697378,
                5224610,
                5223205,
                5821069,
                5740652,
                10026325,
            ),
            (
                "4FB8C44FB3410AA3E90AA41F58D27F57970D98FFCE3261766E0008",
                "2026-10-13T23:08:30",
                8,
                697321,
                697375,
                5224644,
                5223233,
                5821055,
                5740301,
                10026958,
            ),
            # Below is test data provided by Adam N
            (
                "6291E06082EE05DED205C7FE681A7D7FEC287FF4823232B3AC0000",
                "2026-09-08T11:53:16",
                0,
                384722,
                378878,
                6459872,
                6324974,
                6822525,
                8383528,
                8385666,
            ),
        ],
    )
    def test_read_phyre_hex(
        self,
        hex_segment,
        date,
        event,
        ph_temperature_external,
        ph_temperature_internal,
        vrs_external,
        vrs_internal,
        vk,
        ib,
        ik,
    ):
        result = si.read_hex(
            hex_segment,
            si.InstrumentType.pHyre,
            [
                si.SensorID.PHYRE_PH_EXTERNAL,
                si.SensorID.PHYRE_PH_INTERNAL,
                si.SensorID.PHYRE_TEMPERATURE_EXTERNAL,
                si.SensorID.PHYRE_TEMPERATURE_INTERNAL,
                si.SensorID.PH_VK,
                si.SensorID.PH_IB,
                si.SensorID.PH_IK,
                si.SensorID.INSTRUMENT_DATE_TIME,
                si.SensorID.ERROR_FLAG,
            ],
        )

        assert result[str(si.SensorID.PHYRE_PH_EXTERNAL) + "_0"] == vrs_external
        assert result[str(si.SensorID.PHYRE_PH_INTERNAL) + "_0"] == vrs_internal
        assert (
            result[str(si.SensorID.PHYRE_TEMPERATURE_EXTERNAL) + "_0"] == ph_temperature_external
        )
        assert (
            result[str(si.SensorID.PHYRE_TEMPERATURE_INTERNAL) + "_0"] == ph_temperature_internal
        )
        assert result[str(si.SensorID.PH_VK) + "_0"] == vk
        assert result[str(si.SensorID.PH_IB) + "_0"] == ib
        assert result[str(si.SensorID.PH_IK) + "_0"] == ik
        assert result[str(si.SensorID.INSTRUMENT_DATE_TIME) + "_0"] == datetime.fromisoformat(date)
        assert result[str(si.SensorID.ERROR_FLAG) + "_0"] == event


class TestReadpHyreCAT:
    # The following expected values were output directly by a pHyreCAT with the TSS command and
    # outputformat=0. The samples were then uploaded using the uploaddata command in terminal.
    @pytest.mark.parametrize(
        "hex_segment, date, event, ph_temperature_external, ph_temperature_internal, "
        "vrs_external, vrs_internal, vk, ib, ik, "
        "temperature, conductivity, pressure, pressure_temp, ox_phase, ox_temp",
        [
            (
                "0365BB0A45BD08056A0AA728B900181BD54FB88A4FB3320AA3DF0AA42258D3F85799BD98FE5F32616C380008",
                "2026-10-13T22:24:56",
                8,
                697311,
                697378,
                5224586,
                5223218,
                5821432,
                5740989,
                10026591,
                222651,
                673213,
                525674,
                2727,
                2668800,
                1579989,
            ),
            (
                "0365EF0A45BC08056A0AA828B964181BD64FB8B44FB31F0AA3DC0AA42958D3E9579A3598FF5932616C400008",
                "2026-10-13T22:25:04",
                8,
                697308,
                697385,
                5224628,
                5223199,
                5821417,
                5741109,
                10026841,
                222703,
                673212,
                525674,
                2728,
                2668900,
                1579990,
            ),
            (
                "03660F0A45BC0805690AA728B900181BD94FB8884FB3370AA3CE0AA42258D3D457987598FEFA32616C470008",
                "2026-10-13T22:25:11",
                8,
                697294,
                697378,
                5224584,
                5223223,
                5821396,
                5740661,
                10026746,
                222735,
                673212,
                525673,
                2727,
                2668800,
                1579993,
            ),
        ],
    )
    def test_read_phyre_cat_hex(
        self,
        hex_segment,
        date,
        event,
        ph_temperature_external,
        ph_temperature_internal,
        vrs_external,
        vrs_internal,
        vk,
        ib,
        ik,
        temperature,
        conductivity,
        pressure,
        pressure_temp,
        ox_phase,
        ox_temp,
    ):
        result = si.read_hex(
            hex_segment,
            si.InstrumentType.pHyreCAT,
            [
                si.SensorID.TEMPERATURE,
                si.SensorID.CONDUCTIVITY,
                si.SensorID.STRAIN_GAUGE_PRESSURE,
                si.SensorID.STRAIN_GAUGE_PRESSURE_TEMPERATURE,
                si.SensorID.SBE63_DISSOLVED_OXYGEN_PHASE,
                si.SensorID.SBE63_DISSOLVED_OXYGEN_TEMPERATURE,
                si.SensorID.PHYRE_PH_EXTERNAL,
                si.SensorID.PHYRE_PH_INTERNAL,
                si.SensorID.PHYRE_TEMPERATURE_EXTERNAL,
                si.SensorID.PHYRE_TEMPERATURE_INTERNAL,
                si.SensorID.PH_VK,
                si.SensorID.PH_IB,
                si.SensorID.PH_IK,
                si.SensorID.INSTRUMENT_DATE_TIME,
                si.SensorID.ERROR_FLAG,
            ],
        )
        assert result[str(si.SensorID.TEMPERATURE) + "_0"] == temperature
        assert result[str(si.SensorID.CONDUCTIVITY) + "_0"] == conductivity
        assert result[str(si.SensorID.STRAIN_GAUGE_PRESSURE) + "_0"] == pressure
        assert result[str(si.SensorID.STRAIN_GAUGE_PRESSURE_TEMPERATURE) + "_0"] == pressure_temp
        assert result[str(si.SensorID.SBE63_DISSOLVED_OXYGEN_PHASE) + "_0"] == ox_phase
        assert result[str(si.SensorID.SBE63_DISSOLVED_OXYGEN_TEMPERATURE) + "_0"] == ox_temp
        assert result[str(si.SensorID.PHYRE_PH_EXTERNAL) + "_0"] == vrs_external
        assert result[str(si.SensorID.PHYRE_PH_INTERNAL) + "_0"] == vrs_internal
        assert (
            result[str(si.SensorID.PHYRE_TEMPERATURE_EXTERNAL) + "_0"] == ph_temperature_external
        )
        assert (
            result[str(si.SensorID.PHYRE_TEMPERATURE_INTERNAL) + "_0"] == ph_temperature_internal
        )
        assert result[str(si.SensorID.PH_VK) + "_0"] == vk
        assert result[str(si.SensorID.PH_IB) + "_0"] == ib
        assert result[str(si.SensorID.PH_IK) + "_0"] == ik

        assert result[str(si.SensorID.INSTRUMENT_DATE_TIME) + "_0"] == datetime.fromisoformat(date)
        assert result[str(si.SensorID.ERROR_FLAG) + "_0"] == event


class TestReadHexNoDataRows:
    def test_read_hex_file_without_scans(self, tmp_path):
        """A .hex file whose header is not followed by any scans gives an empty dataset."""
        filepath = tmp_path / "header_only.hex"
        filepath.write_text("* Sea-Bird SBE19plus Data File:\n*END*\n")
        raw = si.read_hex_file(filepath, si.InstrumentType.SBE19Plus, [])
        assert len(raw.data_vars) == 0


class TestReadHexIncorrectLength:
    @pytest.mark.parametrize(
        "hex_segment",
        ["00000100002", "0000010000020"],
        ids=["too_short", "too_long"],
    )
    def test_incorrect_hex_length_sets_all_values_to_nan(self, hex_segment):
        sensor_ids = [si.SensorID.TEMPERATURE, si.SensorID.CONDUCTIVITY]
        result = si.read_hex(hex_segment, si.InstrumentType.SBE19Plus, sensor_ids)

        assert set(result) == {str(sensor_id) + "_0" for sensor_id in sensor_ids}
        assert all(pd.isna(value) for value in result.values())
