"""instrument data conversion unit tests."""

from datetime import datetime
from pathlib import Path

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


class TestReadHex19plus:
    filepath = test_data / "19plus_V2_trim.hex"

    @pytest.fixture
    def hex_data(self):
        return si.read_hex_file(
            self.filepath,
            si.InstrumentType.SBE19Plus,
            [
                si.Sensors.Temperature,
                si.Sensors.Conductivity,
                si.Sensors.Pressure,
                si.Sensors.ExtVolt0,
                si.Sensors.ExtVolt1,
                si.Sensors.ExtVolt2,
                si.Sensors.ExtVolt4,
            ],
        )

    def test_read_19plus_hex(self, hex_data):
        assert hex_data["temperature"][0].item() == 322798
        assert hex_data["temperature"][-1].item() == 326317
        assert hex_data["conductivity"][0].item() == 675415 / 256
        assert hex_data["conductivity"][-1].item() == 1590943 / 256
        assert hex_data["pressure"][0].item() == 533539
        assert hex_data["pressure"][-1].item() == 533538
        assert hex_data["temperature compensation"][0].item() == 20625 / 13107
        assert hex_data["temperature compensation"][-1].item() == 20361 / 13107


class TestReadHex16plusDigiquartz:
    filepath = test_data / "16plus_digiquartz.hex"

    @pytest.fixture
    def hex_data(self):
        return si.read_hex_file(
            self.filepath,
            si.InstrumentType.SBE16Plus,
            [
                si.Sensors.Temperature,
                si.Sensors.Conductivity,
                si.Sensors.DigiquartzPressure,
            ],
        )

    def test_read_16plus_hex_digiquartz(self, hex_data):
        assert hex_data["temperature"][0].item() == 407128
        assert hex_data["temperature"][-1].item() == 468996
        assert hex_data["conductivity"][0].item() == 687669 / 256
        assert hex_data["conductivity"][-1].item() == 687669 / 256
        assert hex_data["digiquartz pressure"][0].item() == 8906797 / 256
        assert hex_data["digiquartz pressure"][-1].item() == 8905702 / 256
        assert hex_data["temperature compensation"][0].item() == 1571
        assert hex_data["temperature compensation"][-1].item() == 1566


class TestReadHex39plus:
    filepath = test_data / "SBE39plus.hex"

    @pytest.fixture
    def hex_data(self):
        return si.read_hex_file(
            self.filepath,
            si.InstrumentType.SBE39Plus,
            [
                si.Sensors.Temperature,
                si.Sensors.Conductivity,
                si.Sensors.Pressure,
            ],
        )

    def test_read_39plus_hex(self, hex_data):
        assert hex_data["temperature"][0].item() == 1804057
        assert hex_data["temperature"][-1].item() == 1803910
        assert hex_data["pressure"][0].item() == 108184
        assert hex_data["pressure"][-1].item() == 108132
        assert hex_data["temperature compensation"][0].item() == 445677321 / 13107
        assert hex_data["temperature compensation"][-1].item() == 445664214 / 13107
        assert pd.to_datetime(hex_data["date time"][0].item()) == datetime(2025, 9, 15, 20, 27, 22)
        assert pd.to_datetime(hex_data["date time"][-1].item()) == datetime(2025, 9, 15, 20, 28, 9)


class TestReadSBE911plus:
    filepath = test_data / "911plus_trim.hex"

    @pytest.fixture
    def hex_data(self):
        return si.read_hex_file(
            self.filepath,
            si.InstrumentType.SBE911Plus,
            [
                si.Sensors.Temperature,
                si.Sensors.Conductivity,
                si.Sensors.Pressure,
                si.Sensors.SecondaryTemperature,
                si.Sensors.SecondaryConductivity,
                si.Sensors.ExtVolt0,
                si.Sensors.ExtVolt1,
                si.Sensors.ExtVolt2,
                si.Sensors.ExtVolt3,
                si.Sensors.ExtVolt4,
                si.Sensors.ExtVolt5,
                si.Sensors.ExtVolt6,
                si.Sensors.ExtVolt7,
                si.Sensors.SPAR,
                si.Sensors.nmeaLocation,
                si.Sensors.SystemTime,
            ],
        )

    # Expected values are from SBE Data Processing
    def test_read_911plus_hex(self, hex_data):
        assert round(hex_data["temperature"][0].item(), 3) == 4651.168
        assert round(hex_data["conductivity"][0].item(), 3) == 6638.984
        assert round(hex_data["digiquartz pressure"][0].item(), 3) == 33302.258
        assert round(hex_data["secondary temperature"][0].item(), 3) == 4651.781
        assert round(hex_data["secondary conductivity"][0].item(), 3) == 6230.480
        assert round(hex_data["volt 0"][0].item(), 4) == 3.7338
        assert round(hex_data["volt 1"][0].item(), 4) == 0.2112
        assert round(hex_data["volt 2"][0].item(), 4) == 4.9328
        assert round(hex_data["volt 3"][0].item(), 4) == 2.3370
        assert round(hex_data["volt 4"][0].item(), 4) == 2.4884
        assert round(hex_data["volt 5"][0].item(), 4) == 2.6777
        assert round(hex_data["volt 6"][0].item(), 4) == 0.1502
        assert round(hex_data["volt 7"][0].item(), 4) == 2.7912
        assert round(hex_data["surface par"][0].item(), 5) == 1.43101
        assert hex_data["NMEA Latitude"][0].item() == 34.2757
        assert hex_data["NMEA Longitude"][0].item() == -120.0256
        # Status nibble 0b0011: pump on, no bottom contact (pump and bottom contact match SBE
        # Data Processing; Fathom read the bits most significant first)
        assert int(hex_data["SBE911 pump status"][0].item()) == 1
        assert int(hex_data["SBE911 bottom contact status"][0].item()) == 1
        assert int(hex_data["SBE911 confirm status"][0].item()) == 0
        assert int(hex_data["SBE911 modem status"][0].item()) == 0
        assert int(hex_data["data integrity"][0].item()) == 83
        assert pd.to_datetime(hex_data["system time"][0].item()) == datetime(
            2025, 7, 28, 23, 43, 53
        )


class TestReadSeaFET2:
    filepath = test_data / "SeaFET2.hex"

    @pytest.fixture
    def hex_data(self):
        return si.read_hex_file(self.filepath, si.InstrumentType.SeaFET2, [], False, True)

    def test_read_seaFET2_hex(self, hex_data):
        assert hex_data["vrs external"][0].item() == 5106181
        assert hex_data["vrs internal"][0].item() == 5096276
        assert hex_data["ph temperature"][0].item() == 1002679
        assert hex_data["vk"][0].item() == 5960288
        assert hex_data["ib"][0].item() == 8381912
        assert hex_data["ik"][0].item() == 8384541
        assert hex_data["relative humidity"][0].item() == 17632
        assert hex_data["internal temperature"][0].item() == 26016
        assert hex_data["error flag"][0].item() == 0
        assert pd.to_datetime(hex_data["date time"][0].item()) == datetime(2025, 9, 15, 21, 35, 48)


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
        result = si.read_hex(si.InstrumentType.pHyre, hex_segment)
        assert result[si.HEX_TYPE_DATE_TIME] == datetime.fromisoformat(date)
        assert result[si.HEX_TYPE_ERROR_FLAG] == event
        assert result[si.HEX_TYPE_PH_TEMPERATURE_EXTERNAL] == ph_temperature_external
        assert result[si.HEX_TYPE_PH_TEMPERATURE_INTERNAL] == ph_temperature_internal
        assert result[si.HEX_TYPE_VRS_EXTERNAL] == vrs_external
        assert result[si.HEX_TYPE_VRS_INTERNAL] == vrs_internal
        assert result[si.HEX_TYPE_VK] == vk
        assert result[si.HEX_TYPE_IB] == ib
        assert result[si.HEX_TYPE_IK] == ik


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
        result = si.read_hex(si.InstrumentType.pHyreCAT, hex_segment)
        assert result[si.HEX_TYPE_DATE_TIME] == datetime.fromisoformat(date)
        assert result[si.HEX_TYPE_ERROR_FLAG] == event
        assert result[si.HEX_TYPE_PH_TEMPERATURE_EXTERNAL] == ph_temperature_external
        assert result[si.HEX_TYPE_PH_TEMPERATURE_INTERNAL] == ph_temperature_internal
        assert result[si.HEX_TYPE_VRS_EXTERNAL] == vrs_external
        assert result[si.HEX_TYPE_VRS_INTERNAL] == vrs_internal
        assert result[si.HEX_TYPE_VK] == vk
        assert result[si.HEX_TYPE_IB] == ib
        assert result[si.HEX_TYPE_IK] == ik

        assert result[si.HEX_TYPE_TEMPERATURE] == temperature
        assert result[si.HEX_TYPE_CONDUCTIVITY] == conductivity
        assert result[si.HEX_TYPE_PRESSURE] == pressure
        assert result[si.HEX_TYPE_TEMPERATURE_COMPENSATION] == pressure_temp
        assert result[si.HEX_TYPE_SBE63_PHASE] == ox_phase
        assert result[si.HEX_TYPE_SBE63_TEMPERATURE] == ox_temp

    def test_wrong_instrument_type(self):
        with pytest.raises(ValueError):
            si.read_phyre_data("", si.InstrumentType.SeaFET2)

    def test_wrong_hex_length(self):
        with pytest.raises(ValueError):
            si.read_phyre_data(
                "4FB8A04FB3290AA3EE0AA42258D2A357994E99004B326176670008",
                si.InstrumentType.pHyreCAT,
            )


class TestReadNMEAHex:
    filepath = test_data / "19plus_V2_CTD_nmea_data.hex"

    @pytest.fixture
    def hex_data(self):
        return si.read_hex_file(
            self.filepath,
            si.InstrumentType.SBE19Plus,
            [
                si.Sensors.Temperature,
                si.Sensors.Conductivity,
                si.Sensors.Pressure,
                si.Sensors.ExtVolt0,
                si.Sensors.ExtVolt1,
                si.Sensors.nmeaLatitude,
                si.Sensors.nmeaLongitude,
                si.Sensors.statusAndSign,
                si.Sensors.nmeaTime,
            ],
        )

    def test_read_nmea_hex(self, hex_data):
        assert hex_data["NMEA Latitude"][0].item() == 48.62280000
        assert hex_data["NMEA Latitude"][-1].item() == 48.62298000
        assert hex_data["NMEA Longitude"][0].item() == -123.50002000
        assert hex_data["NMEA Longitude"][-1].item() == -123.49994000
        assert pd.to_datetime(hex_data["NMEA Date Time"][0].item()) == datetime(
            2023, 3, 18, 12, 15, 19
        )
        assert pd.to_datetime(hex_data["NMEA Date Time"][-1].item()) == datetime(
            2023, 3, 18, 12, 36, 2
        )


class TestReadHexNoDataRows:
    def test_read_hex_file_without_scans(self, tmp_path):
        """A .hex file whose header is not followed by any scans gives an empty dataset."""
        filepath = tmp_path / "header_only.hex"
        filepath.write_text("* Sea-Bird SBE19plus Data File:\n*END*\n")
        raw = si.read_hex_file(filepath, si.InstrumentType.SBE19Plus, [], False, True)
        assert len(raw.data_vars) == 0


class TestSBE911PlusStatusBits:
    # 5 frequencies, 4 voltage words, pressure sensor temperature, status, modulo count
    enabled_sensors = (
        si.Sensors.SecondaryTemperature,
        si.Sensors.SecondaryConductivity,
        *(getattr(si.Sensors, f"ExtVolt{n}") for n in range(8)),
    )

    def scan(self, nibble: int) -> str:
        return "0" * 30 + "0" * 24 + "ABC" + f"{nibble:X}" + "6C"

    @pytest.mark.parametrize("nibble", range(16))
    def test_bits_read_least_significant_first(self, nibble):
        result = si.read_sbe911plus_data(self.scan(nibble), self.enabled_sensors)
        assert result["SBE911 pump status"] == nibble & 1
        assert result["SBE911 bottom contact status"] == (nibble >> 1) & 1
        assert result["SBE911 confirm status"] == (nibble >> 2) & 1
        assert result["SBE911 modem status"] == (nibble >> 3) & 1
        assert result["temperature compensation"] == 0xABC
        assert result["data integrity"] == 0x6C

    @pytest.mark.parametrize(
        "nibble, pump, confirm",
        [
            (0b0010, 0, 0),  # deck test: pump off
            (0b0011, 1, 0),  # in the water: pump on
            (0b0111, 1, 1),  # bottle firing
        ],
    )
    def test_status_seen_in_real_casts(self, nibble, pump, confirm):
        result = si.read_sbe911plus_data(self.scan(nibble), self.enabled_sensors)
        assert result["SBE911 pump status"] == pump
        assert result["SBE911 confirm status"] == confirm
        # no bottom contact, modem carrier detected
        assert result["SBE911 bottom contact status"] == 1
        assert result["SBE911 modem status"] == 0
