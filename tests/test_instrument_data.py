"""instrument data conversion unit tests."""

# Native imports
from datetime import datetime
from pathlib import Path

# Third-party imports
import pandas as pd
import pytest

# Sea-Bird imports

# Internal imports
import seabirdscientific.instrument_data as si

test_data = Path("./tests/resources/test-data")


class TestCnvToInstrumentData:
    def test_read_cnv_file_pass(self):
        file_path = test_data / "SBE37SM.cnv"
        dataset = si.read_cnv_file(file_path)

        expected_variables = set(
            [
                "tv290C",
                "cond0S/m",
                "prdM",
                "prdE",
                "depSM",
                "depSF",
                "sal00",
                "flag",
            ]
        )
        assert expected_variables == set(dataset.data_vars)
        assert dataset.attrs["sample_interval"] == 120
        assert dataset["tv290C"].name == "tv290C"
        assert dataset["tv290C"].attrs["long_name"] == "Temperature"
        assert dataset["tv290C"].attrs["units"] == "ITS-90, deg C"
        assert dataset["tv290C"][0] == 20.8989
        assert dataset["tv290C"][-1] == -45.7713

    def test_duplicate_labels(self):
        """
        Verify that duplicate labels are handled correctly.
        This .cnv file is from the Australian Institute of Marine Science
        """
        file_path = test_data / "WQR086CFACLWDB.cnv"
        dataset = si.read_cnv_file(file_path)

        expected_variables = set(
            [
                "prDM",
                "t090C",
                "c0S/m",
                "scan",
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
            ]
        )

        assert expected_variables == set(dataset.data_vars)

        # MeasurementSeries should have the original label
        assert dataset["depSM_1"].attrs["sbs_name"] == "depSM"

        first_values = [dataset[key].sel(sample=0).item() for key in list(dataset.data_vars)]
        first_data_line = "1.006 28.0831 5.811991 5248 0.998 218.616 6.9285e+01 3.9379e+02 11.8057 0.8223 2.7316e+02 28.0829 5.811250 1.000 36.2666 31 0.0000e+00"
        expected_first_values = [float(v) for v in first_data_line.split()]
        assert first_values == expected_first_values


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
        # Does not match SBE Processing, matches Fathom
        assert int(hex_data["SBE911 pump status"][0].item()) == 0
        # Does not match SBE Processing, matches Fathom
        assert int(hex_data["SBE911 bottom contact status"][0].item()) == 0
        assert int(hex_data["SBE911 confirm status"][0].item()) == 1
        assert int(hex_data["SBE911 modem status"][0].item()) == 1
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
