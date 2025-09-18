"""instrument data conversion unit tests."""

# Native imports
from datetime import datetime
from pathlib import Path

# Third-party imports

# Sea-Bird imports

# Internal imports
import seabirdscientific.instrument_data as id

test_data = Path("./resources/test-data")


class TestCnvToInstrumentData:
    measurement_series = id.MeasurementSeries(
        label="test",
        description="",
        units="",
        start_time=datetime(2000, 1, 1, 1, 1, 1),
        values=[10, 20, 30],
    )

    def test_create_measurement_series_pass(self):
        assert self.measurement_series.label == "test"

    def test_create_instrument_data_pass(self):
        data = id.InstrumentData(
            measurements={"test": self.measurement_series},
            interval_s=None,
            latitude=0.0,
            start_time=None,
            sample_count=None,
        )
        assert True  # fails before assert if data fails to instantiate

    def test_cnv_to_instrument_data_pass(self):
        file_path = test_data / "SBE37SM.cnv"
        data = id.cnv_to_instrument_data(file_path)

        expected_labels = [
            "tv290C",
            "cond0S/m",
            "prdM",
            "prdE",
            "depSM",
            "depSF",
            "sal00",
            "flag",
        ]
        for expected, result in zip(expected_labels, data.measurements):
            assert expected == result

        assert data.interval_s == 120

        assert data.measurements["tv290C"].label == "tv290C"
        assert data.measurements["tv290C"].description == "Temperature"
        assert data.measurements["tv290C"].units == "ITS-90, deg C"
        assert data.measurements["tv290C"].start_time == datetime(2014, 2, 19, 13, 30, 1)
        assert data.measurements["tv290C"].values[0] == 20.8989
        assert data.measurements["tv290C"].values[-1] == -45.7713

    def test_duplicate_labels(self):
        """
        Verify that duplicate labels are handled correctly.
        This .cnv file is from the Australian Institute of Marine Science
        """
        file_path = test_data / "WQR086CFACLWDB.cnv"
        data = id.cnv_to_instrument_data(file_path)

        expected_labels = [
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
            "depSM_1",  # this duplicate "depSM" was a problem, increments to "depSM1"
            "sal00",
            "nbin",
            "flag",
        ]

        assert list(data.measurements.keys()) == expected_labels

        # MeasurementSeries should have the original label
        assert data.measurements["depSM_1"].label == "depSM"

        first_values = [m.values[0] for m in data.measurements.values()]
        first_data_line = "1.006 28.0831 5.811991 5248 0.998 218.616 6.9285e+01 3.9379e+02 11.8057 0.8223 2.7316e+02 28.0829 5.811250 1.000 36.2666 31 0.0000e+00"
        expected_first_values = [float(v) for v in first_data_line.split()]
        assert first_values == expected_first_values


class TestReadHex19plus:
    filepath = test_data / "19plus_V2" / "19plus_V2.hex"
    raw = id.read_hex_file(
        filepath,
        id.InstrumentType.SBE19Plus,
        [
            id.Sensors.Temperature,
            id.Sensors.Conductivity,
            id.Sensors.Pressure,
            id.Sensors.ExtVolt0,
            id.Sensors.ExtVolt1,
            id.Sensors.ExtVolt2,
            id.Sensors.ExtVolt4,
        ],
    )

    def test_read_19plus_hex(self):
        assert self.raw["temperature"].iloc[0] == 322798
        assert self.raw["temperature"].iloc[-1] == 326317
        assert self.raw["conductivity"].iloc[0] == 675415 / 256
        assert self.raw["conductivity"].iloc[-1] == 1590943 / 256
        assert self.raw["pressure"].iloc[0] == 533539
        assert self.raw["pressure"].iloc[-1] == 533538
        assert self.raw["temperature compensation"].iloc[0] == 20625 / 13107
        assert self.raw["temperature compensation"].iloc[-1] == 20361 / 13107


class TestReadHex39plus:
    filepath = test_data / "SBE39plus" / "SBE39plus.hex"
    raw = id.read_hex_file(
        filepath,
        id.InstrumentType.SBE39Plus,
        [
            id.Sensors.Temperature,
            id.Sensors.Conductivity,
            id.Sensors.Pressure,
        ],
    )

    def test_read_39plus_hex(self):
        assert self.raw["temperature"].iloc[0] == 1804057
        assert self.raw["temperature"].iloc[-1] == 1803910
        assert self.raw["pressure"].iloc[0] == 108184
        assert self.raw["pressure"].iloc[-1] == 108132
        assert self.raw["temperature compensation"].iloc[0] == 445677321 / 13107
        assert self.raw["temperature compensation"].iloc[-1] == 445664214 / 13107
        assert self.raw["date time"].iloc[0] == datetime(2025, 9, 15, 20, 27, 22)
        assert self.raw["date time"].iloc[-1] == datetime(2025, 9, 15, 20, 28, 9)


class TestReadSBE911plus:
    filepath = test_data / "SBE911plus/2507054.hex"
    raw = id.read_hex_file(
        filepath,
        id.InstrumentType.SBE911Plus,
        [
            id.Sensors.Temperature,
            id.Sensors.Conductivity,
            id.Sensors.Pressure,
            id.Sensors.SecondaryTemperature,
            id.Sensors.SecondaryConductivity,
            id.Sensors.ExtVolt0,
            id.Sensors.ExtVolt1,
            id.Sensors.ExtVolt2,
            id.Sensors.ExtVolt3,
            id.Sensors.ExtVolt4,
            id.Sensors.ExtVolt5,
            id.Sensors.ExtVolt6,
            id.Sensors.ExtVolt7,
            id.Sensors.SPAR,
            id.Sensors.nmeaLocation,
            id.Sensors.SystemTime,
        ],
    )

    # Expected values are from SBE Data Processing
    def test_read_911plus_hex(self):
        print(self.raw.iloc[0])
        assert round(self.raw["temperature"].iloc[0], 3) == 4651.168
        assert round(self.raw["conductivity"].iloc[0], 3) == 6638.984
        assert round(self.raw["digiquartz pressure"].iloc[0], 3) == 33302.258
        assert round(self.raw["secondary temperature"].iloc[0], 3) == 4651.781
        assert round(self.raw["secondary conductivity"].iloc[0], 3) == 6230.480
        assert round(self.raw["volt 0"].iloc[0], 4) == 3.7338
        assert round(self.raw["volt 1"].iloc[0], 4) == 0.2112
        assert round(self.raw["volt 2"].iloc[0], 4) == 4.9328
        assert round(self.raw["volt 3"].iloc[0], 4) == 2.3370
        assert round(self.raw["volt 4"].iloc[0], 4) == 2.4884
        assert round(self.raw["volt 5"].iloc[0], 4) == 2.6777
        assert round(self.raw["volt 6"].iloc[0], 4) == 0.1502
        assert round(self.raw["volt 7"].iloc[0], 3) == 2.7912
        assert round(self.raw["surface par"].iloc[0], 5) == 1.4310
        assert self.raw["NMEA Latitude"].iloc[0] == 34.2757
        assert self.raw["NMEA Longitude"].iloc[0] == -120.0256
        assert (
            int(self.raw["SBE911 pump status"].iloc[0]) == 0
        )  # Does not match SBE Processing, matches Fathom
        assert (
            int(self.raw["SBE911 bottom contact status"].iloc[0]) == 0
        )  # Does not match SBE Processing, matches Fathom
        assert int(self.raw["SBE911 confirm status"].iloc[0]) == 1
        assert int(self.raw["SBE911 modem status"].iloc[0]) == 1
        assert int(self.raw["data integrity"].iloc[0]) == 83
        assert self.raw["system time"].iloc[0] == datetime(2025, 7, 28, 23, 43, 53)


class TestReadSeaFET2:
    filepath = test_data / "SeaFET_V2/SeaFET2.hex"
    raw = id.read_hex_file(filepath, id.InstrumentType.SeaFET2, [], False, True)

    def test_read_seaFET2_hex(self):
        assert self.raw["vrs external"].iloc[0] == 5106181
        assert self.raw["vrs internal"].iloc[0] == 5096276
        assert self.raw["ph temperature"].iloc[0] == 1002679
        assert self.raw["vk"].iloc[0] == 5960288
        assert self.raw["ib"].iloc[0] == 8381912
        assert self.raw["ik"].iloc[0] == 8384541
        assert self.raw["relative humidity"].iloc[0] == 17632
        assert self.raw["internal temperature"].iloc[0] == 26016
        assert self.raw["error flag"].iloc[0] == 0
        assert self.raw["date time"].iloc[0] == datetime(2025, 9, 15, 21, 35, 48)


class TestReadNMEAHex:
    filepath = test_data / "19plus_V2_CTD_nmea_data.hex"
    raw = id.read_hex_file(
        filepath,
        id.InstrumentType.SBE19Plus,
        [
            id.Sensors.Temperature,
            id.Sensors.Conductivity,
            id.Sensors.Pressure,
            id.Sensors.ExtVolt0,
            id.Sensors.ExtVolt1,
            id.Sensors.nmeaLatitude,
            id.Sensors.nmeaLongitude,
            id.Sensors.statusAndSign,
            id.Sensors.nmeaTime,
        ],
    )

    def test_read_nmea_hex(self):
        assert self.raw["NMEA Latitude"].iloc[0] == 48.62280000
        assert self.raw["NMEA Latitude"].iloc[-1] == 48.62298000
        assert self.raw["NMEA Longitude"].iloc[0] == -123.50002000
        assert self.raw["NMEA Longitude"].iloc[-1] == -123.49994000
        assert self.raw["NMEA Date Time"].iloc[0] == datetime(2023, 3, 18, 12, 15, 19)
        assert self.raw["NMEA Date Time"].iloc[-1] == datetime(2023, 3, 18, 12, 36, 2)
