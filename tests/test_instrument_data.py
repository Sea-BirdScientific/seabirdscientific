"""instrument data conversion unit tests.

"""

# Native imports
from datetime import datetime
from pathlib import Path

# Third-party imports

# Sea-Bird imports

# Internal imports
import seabirdscientific.instrument_data as id

test_data = Path("./tests/resources/test-data")


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


class TestReadHex:
    filepath = test_data / "19plus_V2_CTD-processing_example.hex"
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

    def test_read_hex(self):
        assert self.raw["temperature"].iloc[0] == 322798
        assert self.raw["temperature"].iloc[-1] == 326317
        assert self.raw["conductivity"].iloc[0] == 675415 / 256
        assert self.raw["conductivity"].iloc[-1] == 1590943 / 256
        assert self.raw["pressure"].iloc[0] == 533539
        assert self.raw["pressure"].iloc[-1] == 533538
        assert self.raw["temperature compensation"].iloc[0] == 20625 / 13107
        assert self.raw["temperature compensation"].iloc[-1] == 20361 / 13107


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
            id.Sensors.nmeaTime
        ],
    )

    def test_read_nmea_hex(self):
        assert self.raw["NMEA Latitude"].iloc[0] == 48.62280000
        assert self.raw["NMEA Latitude"].iloc[-1] == 48.62298000
        assert self.raw["NMEA Longitude"].iloc[0] == -123.50002000
        assert self.raw["NMEA Longitude"].iloc[-1] == -123.49994000
        assert self.raw["NMEA Date Time"].iloc[0] == datetime.fromisoformat("2023-03-18 05:15:19")
        assert self.raw["NMEA Date Time"].iloc[-1] == datetime.fromisoformat("2023-03-18 05:36:02")
