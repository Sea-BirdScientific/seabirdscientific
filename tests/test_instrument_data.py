"""instrument data conversion unit tests.

"""

# Native imports
from datetime import datetime

# Third-party imports

# Sea-Bird imports

# Internal imports
import seabirdscientific.instrument_data as id


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
        file_path = "./tests/resources/test-data/SBE37SM.cnv"
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


class TestReadHex:
    filepath = "./documentation/example_data/19plus_V2_CTD-processing_example.hex"
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
