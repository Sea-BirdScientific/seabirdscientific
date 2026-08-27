"""Data conversion unit tests."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import seabirdscientific.constants as const
import seabirdscientific.conversion as dc
import seabirdscientific.eos80_conversion as eos80
import seabirdscientific.instrument_data as id
import test_coefficients as tc

test_data = Path("./tests/resources/test-data")


class TestConvertTemperature:
    # test temperature raw values
    test_temp_vals = np.array([322798, 322808, 322827, 322838])
    test_temp_no_mv_r = np.array([2864.51635696, 2864.61073548, 2864.79005751, 2864.89387722])

    def test_convert_temperature_90C(self, request):
        expected = [20.4459, 20.4451, 20.4436, 20.4427]
        result = dc.convert_temperature(
            self.test_temp_vals,
            tc.temperature_coefs_sn6130,
            "ITS90",
            "C",
            True,
        )
        request.node.return_value = result.tolist()
        assert np.allclose(expected, result, rtol=0, atol=1e-4)

    def test_convert_temperature_90F(self, request):
        expected = [68.8027, 68.8012, 68.7984, 68.7968]
        result = dc.convert_temperature(
            self.test_temp_vals,
            tc.temperature_coefs_sn6130,
            "ITS90",
            "F",
            True,
        )
        request.node.return_value = result.tolist()
        assert np.allclose(expected, result, rtol=0, atol=1e-4)

    def test_convert_temperature_68C(self, request):
        expected = [20.4508, 20.4500, 20.4485, 20.4476]
        result = dc.convert_temperature(
            self.test_temp_vals,
            tc.temperature_coefs_sn6130,
            "IPTS68",
            "C",
            True,
        )
        request.node.return_value = result.tolist()
        assert np.allclose(expected, result, rtol=0, atol=1e-4)

    def test_convert_temperature_68F(self, request):
        expected = [68.8115, 68.8100, 68.8073, 68.8057]
        result = dc.convert_temperature(
            self.test_temp_vals,
            tc.temperature_coefs_sn6130,
            "IPTS68",
            "F",
            True,
        )
        request.node.return_value = result.tolist()
        assert np.allclose(expected, result, rtol=0, atol=1e-4)

    def test_convert_temperature_90C_no_mv_r(self, request):
        expected = [20.4459, 20.4451, 20.4436, 20.4427]
        result = dc.convert_temperature(
            self.test_temp_no_mv_r,
            tc.temperature_coefs_sn6130,
            "ITS90",
            "C",
            False,
        )
        request.node.return_value = result.tolist()
        assert np.allclose(expected, result, rtol=0, atol=1e-4)


class TestConvertTemperatureFrequency:
    test_temp_freq_vals = np.array(
        [4651.168, 4650.961, 4650.219, 4649.523, 4649.418, 4649.879, 3771.953, 4084.508]
    )

    def test_convert_temperature_frequency_90C(self, request):
        expected = [16.6799, 16.6777, 16.6698, 16.6624, 16.6613, 16.6662, 6.6247, 10.3732]
        result = dc.convert_temperature_frequency(
            self.test_temp_freq_vals,
            tc.temperature_frequency_coefs_sn5102,
            "ITS90",
            "C",
        )
        request.node.return_value = result.tolist()
        assert np.allclose(expected, result, rtol=0, atol=1e-4)

    def test_convert_temperature_frequency_90F(self, request):
        expected = [62.0238, 62.0199, 62.0056, 61.9923, 61.9903, 61.9991, 43.9244, 50.6718]
        result = dc.convert_temperature_frequency(
            self.test_temp_freq_vals,
            tc.temperature_frequency_coefs_sn5102,
            "ITS90",
            "F",
        )
        request.node.return_value = result.tolist()
        assert np.allclose(expected, result, rtol=0, atol=1e-4)

    def test_convert_temperature_frequency_68C(self, request):
        expected = [16.6839, 16.6817, 16.6738, 16.6664, 16.6653, 16.6702, 6.6263, 10.3757]
        result = dc.convert_temperature_frequency(
            self.test_temp_freq_vals,
            tc.temperature_frequency_coefs_sn5102,
            "IPTS68",
            "C",
        )
        request.node.return_value = result.tolist()
        assert np.allclose(expected, result, rtol=0, atol=1e-4)

    def test_convert_temperature_frequency_68F(self, request):
        expected = [62.0310, 62.0271, 62.0128, 61.9995, 61.9975, 62.0063, 43.9273, 50.6763]
        result = dc.convert_temperature_frequency(
            self.test_temp_freq_vals,
            tc.temperature_frequency_coefs_sn5102,
            "IPTS68",
            "F",
        )
        request.node.return_value = result.tolist()
        assert np.allclose(expected, result, rtol=0, atol=1e-4)


class TestConvertPressure:
    # test pressure raw values
    test_pressure_vals = np.array([533539, 533538, 533540, 533537])

    # test pressure temperature compensation raw values
    test_compensation_vals = np.array([20625, 20626, 20626, 20626]) / 13107

    def test_convert_pressure_psia(self, request):
        expected = [14.547, 14.546, 14.549, 14.544]
        result = dc.convert_pressure(
            self.test_pressure_vals,
            self.test_compensation_vals,
            tc.pressure_coefs_sn6130,
            "psia",
        )
        request.node.return_value = result.tolist()
        assert np.allclose(expected, result, rtol=0, atol=1e-3)

    def test_convert_pressure_dbar(self, request):
        expected = [-0.105, -0.106, -0.104, -0.107]
        result = dc.convert_pressure(
            self.test_pressure_vals,
            self.test_compensation_vals,
            tc.pressure_coefs_sn6130,
            "dbar",
        )
        request.node.return_value = result.tolist()
        assert np.allclose(expected, result, rtol=0, atol=1e-3)


class TestDigiquartzPressure:
    # test pressure raw values
    test_pressure_vals = np.array(
        [33302.258, 33302.250, 33302.230, 33302.285, 33302.250, 33302.258]
    )

    # test pressure temperature compensation raw values
    test_compensation_vals = np.array([2498.0, 2499.0, 2498.0, 2499.0, 2499.0, 2498.0])

    def test_convert_pressure_digiquartz_dbar(self, request):
        expected = [2.099, 2.085, 2.051, 2.146, 2.085, 2.099]
        result = dc.convert_pressure_digiquartz(
            self.test_pressure_vals,
            self.test_compensation_vals,
            tc.pressure_digiquartz_coefs_sn5102,
            "dbar",
            0.25,
        )
        # apply slope and offset values
        result = result * 1.00001297 - 2.46118
        request.node.return_value = result.tolist()
        assert np.allclose(expected, result, rtol=0, atol=1e-3)


class TestConductivity19plus:
    cnv_path = test_data / "19plus_V2_CTD-processing_example.cnv"
    hex_path = test_data / "19plus_V2.hex"

    def test_convert_conductivity(self):
        # expected_data = id.read_cnv_file(self.cnv_path)

        raw = id.read_hex_file(
            self.hex_path,
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

        # test_cond_data = np.array([675415, 675404, 675405, 675398, 675391, 675379])

        # test_temp_data = np.array([20.4459, 20.4451, 20.4436, 20.4427, 20.4413, 20.4401])
        temperature = dc.convert_temperature(
            raw["temperature"][0:6].values,
            tc.temperature_coefs_sn6130,
            "IPTS68",
            "C",
            True,
        )

        # test_press_data = np.array([-0.105, -0.106, -0.104, -0.107, -0.105, -0.104])
        pressure = dc.convert_pressure(
            raw["pressure"][0:6].values,
            raw["temperature compensation"][0:6].values,
            tc.pressure_coefs_sn6130,
            "psia",
        )
        expected = [0.008453, 0.008420, 0.008423, 0.008402, 0.008380, 0.008344]
        result = dc.convert_conductivity(
            raw["conductivity"][0:6].values,
            temperature,
            pressure,
            tc.conductivity_coefs_sn6130,
            id.InstrumentType.SBE19Plus,
        )

        assert np.allclose(expected, result, rtol=0, atol=1e-6)


class TestConductivity37SM:
    cnv_path = test_data / "SBE37SM-RS232_03716125_2017_11_16.cnv"
    hex_path = test_data / "SBE37SM-RS232_03716125_2017_11_16.hex"

    def test_convert_conductivity(self):
        # expected_data = id.read_cnv_file(self.cnv_path)

        raw = id.read_hex_file(
            self.hex_path,
            id.InstrumentType.SBE37SM,
            [
                id.Sensors.Temperature,
                id.Sensors.Conductivity,
                id.Sensors.Pressure,
            ],
        )

        temperature = dc.convert_temperature(
            raw["temperature"][0:6].values,
            tc.temperature_coefs_sn6130,
            "IPTS68",
            "C",
            True,
        )

        pressure = dc.convert_pressure(
            raw["pressure"][0:6].values,
            raw["temperature compensation"][0:6].values,
            tc.pressure_coefs_sn16125,
            "psia",
        )
        expected = [2.711842, 2.715786, 2.715857, 2.715846, 2.715846, 2.715857]
        result = dc.convert_conductivity(
            raw["conductivity"][0:6].values,
            temperature,
            pressure,
            tc.conductivity_coefs_sn16125,
            id.InstrumentType.SBE37SM,
        )

        assert np.allclose(expected, result, rtol=0, atol=1e-4)


class TestDeriveDensity:
    data_path = test_data / "SBE37SM-derived.asc"

    @pytest.fixture
    def data(self):
        return pd.read_csv(self.data_path)

    @pytest.mark.parametrize(
        "reference_pressure, expected_column",
        [
            (0.0, "gsw_sigma0A0"),
            (1000.0, "gsw_sigma1A0"),
            (2000.0, "gsw_sigma2A0"),
            (3000.0, "gsw_sigma3A0"),
            (4000.0, "gsw_sigma4A0"),
        ],
    )
    def test_derive_potential_density_from_t_s_p_pass(
        self, data, request, reference_pressure, expected_column
    ):
        temperature = data["t090C"].values
        salinity = data["sal00"].values
        pressure = data["prM"].values
        expected = data[expected_column].values
        result = dc.potential_density_from_t_s_p(
            temperature,
            salinity,
            pressure,
            reference_pressure=reference_pressure,
        )
        request.node.return_value = result.tolist()
        # TODO: improve passing condition for all instances of allclose
        assert np.allclose(result, expected, atol=0.1)

    def test_derive_density_from_t_s_p_pass(self, data, request):
        temperature = data["t090C"].values
        salinity = data["sal00"].values
        pressure = data["prM"].values
        expected = data["gsw_densityA0"].values
        result = dc.density_from_t_s_p(temperature, salinity, pressure)

        request.node.return_value = result.tolist()
        assert np.allclose(result, expected, atol=0.1)

    @pytest.mark.parametrize(
        "reference_pressure, expected_column",
        [
            (0.0, "gsw_sigma0A0"),
            (1000.0, "gsw_sigma1A0"),
            (2000.0, "gsw_sigma2A0"),
            (3000.0, "gsw_sigma3A0"),
            (4000.0, "gsw_sigma4A0"),
        ],
    )
    def test_derive_potential_density_from_t_c_p_pass(
        self, data, request, reference_pressure, expected_column
    ):
        temperature = data["t090C"].values
        conductivity = data["c0S/m"].values * 10.0
        pressure = data["prM"].values
        expected = data[expected_column].values
        result = dc.potential_density_from_t_c_p(
            temperature,
            conductivity,
            pressure,
            reference_pressure=reference_pressure,
        )
        request.node.return_value = result.tolist()
        assert np.allclose(result, expected, atol=0.1)

    def test_derive_density_from_t_c_p_pass(self, data, request):
        temperature = data["t090C"].values
        conductivity = data["c0S/m"].values * 10.0
        pressure = data["prM"].values
        expected = data["gsw_densityA0"].values
        result = dc.density_from_t_c_p(temperature, conductivity, pressure)

        request.node.return_value = result.tolist()
        assert np.allclose(result, expected, atol=0.1)


class TestDepthFromPressure:
    data_path = test_data / "SBE37SM.asc"

    @pytest.fixture
    def data(self):
        return pd.read_csv(self.data_path).loc[6:10]

    @pytest.mark.parametrize(
        "pressure, pressure_units, depth, depth_units",
        [
            ("prdM", "dbar", "depSM", "m"),
            ("prdM", "dbar", "depSF", "ft"),
            ("prdE", "psi", "depSM", "m"),
            ("prdE", "psi", "depSF", "ft"),
        ],
    )
    def test_depth_from_pressure_pass(
        self, data, request, pressure, pressure_units, depth, depth_units
    ):
        expected_depth = data[depth].values
        pressure = data[pressure].values
        result_depth = dc.depth_from_pressure(pressure, 0, depth_units, pressure_units)
        request.node.return_value = result_depth.tolist()
        assert np.allclose(expected_depth, result_depth, atol=0.002)

class TestDeriveSoundVelocity:
    cnv_path = test_data / "SBE19plus_derive_testing.cnv"

    @pytest.fixture
    def source_data(self):
        return id.read_cnv_file(self.cnv_path, "seasoft")

    def test_derive_sound_velocity_c(self, source_data):
        result = eos80.derive_sound_velocity_c(
            source_data["sal00"].values, source_data["tv290C"].values, source_data["prdM"].values
        )
        # TODO: Update expected column name once verified in CNV file
        expected = source_data["svCM"].values

        assert np.allclose(result, expected, rtol=0, atol=1e-1)

    def test_derive_sound_velocity_d(self, source_data):
        result = eos80.derive_sound_velocity_d(
            source_data["sal00"].values, source_data["tv290C"].values, source_data["prdM"].values
        )
        # TODO: Update expected column name once verified in CNV file
        expected = source_data["svDM"].values

        assert np.allclose(result, expected, rtol=0, atol=1e-1)

    def test_derive_sound_velocity_w(self, source_data):
        result = eos80.derive_sound_velocity_w(
            source_data["sal00"].values, source_data["tv290C"].values, source_data["prdM"].values
        )
        # TODO: Update expected column name once verified in CNV file
        expected = source_data["svWM"].values

        assert np.allclose(result, expected, rtol=0, atol=1e-1)


class TestConvertSBE43Oxygen:
    # Some tests need to be run on complete datasets
    cnv_path = test_data / "SBE19plus_01906398_2019_07_15_0033-seasoft-convert-o2-full.cnv"

    @pytest.fixture
    def source_data(self):
        return id.read_cnv_file(self.cnv_path, "seasoft")

    def test_convert_sbe43_oxygen(self, request):
        # From O3287.pdf in the shared calibration folder
        raw_oxygen = np.array([0.725, 0.756, 0.803, 0.874, 0.925, 0.96, 1.332, 1.435, 1.595, 1.81])
        pressure = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        temperature = np.array([2, 6, 12, 20, 26, 30, 2, 6, 12, 20])
        salinity = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        expected = np.array([1.11, 1.12, 1.13, 1.16, 1.18, 1.18, 3.9, 3.9, 3.92, 3.95])

        result = dc._convert_sbe43_oxygen(
            raw_oxygen,
            temperature,
            pressure,
            salinity,
            tc.oxygen_43_coefs_sn3287,
            0,
        )
        request.node.return_value = result.tolist()
        assert np.allclose(expected, result, rtol=0, atol=1e-2)

    def test_convert_sbe43_oxygen_from_hex(self, request):
        # fmt: off
        # From SBE19plus_01906398_2019_07_15_0033.hex
        raw_oxygen = np.array(
            [2.5575, 2.5586, 2.5606, 2.5627, 2.5638, 2.5637, 2.5635, 2.5629, 2.5621, 2.5618]
        )
        pressure = np.array(
            [-0.012, -0.012, -0.012, -0.012, -0.012, -0.012, -0.012, -0.011, 0.107, 0.351]
        )
        temperature = np.array(
            [25.3427, 25.3408, 25.3387, 25.3363, 25.3341, 25.3326, 25.3316, 25.3302, 25.3377, 25.5433]
        )
        salinity = np.array(
            [0.4373, 0.5592, 0.5865, 0.5095, 0.4621, 0.4119, 0.3936, 0.3463, 4.9297, 6.5098]
        )
        expected = np.array(
            [4.4728, 4.4722, 4.4762, 4.4828, 4.4867, 4.4879, 4.488, 4.488, 4.3707, 4.3148]
        )
        # fmt: on
        result = dc._convert_sbe43_oxygen(
            raw_oxygen,
            temperature,
            pressure,
            salinity,
            tc.oxygen_43_coefs_sn1686,
            0,
        )
        request.node.return_value = result.tolist()
        assert np.allclose(expected, result, rtol=0, atol=1e-3)

    def test_convert_sbe43_oxygen_from_hex_with_hysteresis(self, request):
        # fmt: off
        # From SBE19plus_01906398_2019_07_15_0033.hex
        # TODO: hysteresis correction only has a real impact on deep data, will need some to better validate this
        raw_oxygen = np.array(
            [2.5575, 2.5586, 2.5606, 2.5627, 2.5638, 2.5637, 2.5635, 2.5629, 2.5621, 2.5618]
        )
        pressure = np.array(
            [-0.012, -0.012, -0.012, -0.012, -0.012, -0.012, -0.012, -0.011, 0.107, 0.351]
        )
        temperature = np.array(
            [25.3427, 25.3408, 25.3387, 25.3363, 25.3341, 25.3326, 25.3316, 25.3302, 25.3377, 25.5433]
        )
        salinity = np.array(
            [0.4373, 0.5592, 0.5865, 0.5095, 0.4621, 0.4119, 0.3936, 0.3463, 4.9297, 6.5098]
        )
        expected = np.array(
            [4.4728, 4.4722, 4.4762, 4.4828, 4.4867, 4.4879, 4.488, 4.488, 4.3707, 4.3148]
        )
        # fmt: on
        result = dc.convert_sbe43_oxygen(
            raw_oxygen,
            temperature,
            pressure,
            salinity,
            tc.oxygen_43_coefs_sn1686,
            False,
            True,
            1,
            0.25,
        )

        request.node.return_value = result.tolist()
        assert np.allclose(result, expected, rtol=0, atol=1e-3)

    def test_convert_sbe43_oxygen_from_hex_with_tau_correction(self, request):
        # fmt: off
        # From SBE19plus_01906398_2019_07_15_0033.hex
        raw_oxygen = np.asarray(
            [2.5575, 2.5586, 2.5606, 2.5627, 2.5638, 2.5637, 2.5635, 2.5629, 2.5621, 2.5618]
        )
        pressure = np.asarray(
            [-0.012, -0.012, -0.012, -0.012, -0.012, -0.012, -0.012, -0.011, 0.107, 0.351]
        )
        temperature = np.asarray(
            [25.3427, 25.3408, 25.3387, 25.3363, 25.3341, 25.3326, 25.3316, 25.3302, 25.3377, 25.5433]
        )
        salinity = np.asarray(
            [0.4373, 0.5592, 0.5865, 0.5095, 0.4621, 0.4119, 0.3936, 0.3463, 4.9297, 6.5098]
        )
        expected = [4.4729, 4.4723, 4.4884, 4.4927, 4.4916, 4.4879, 4.4849, 4.4841, 4.3707, 4.3148]
        # fmt: on
        result = dc.convert_sbe43_oxygen(
            raw_oxygen,
            temperature,
            pressure,
            salinity,
            tc.oxygen_43_coefs_sn1686,
            True,
            False,
            1,
            0.25,
        )

        request.node.return_value = result.tolist()
        assert np.allclose(result, expected, rtol=0, atol=1e-4)

    def test_convert_to_mg_per_l(self, request):
        oxMlPerL = np.array(
            [4.4728, 4.4722, 4.4762, 4.4828, 4.4867, 4.4879, 4.488, 4.488, 4.3707, 4.3148]
        )
        expected = [6.3921, 6.3913, 6.3969, 6.4064, 6.4119, 6.4137, 6.4138, 6.4138, 6.2461, 6.1663]
        result = dc.convert_oxygen_to_mg_per_l(oxMlPerL)
        request.node.return_value = result.tolist()
        assert np.allclose(expected, result, rtol=0, atol=1e-3)

    def test_convert_to_umol_per_kg(self, request):
        # fmt: off
        oxMlPerL = np.array(
            [4.4728, 4.4722, 4.4762, 4.4828, 4.4867, 4.4879, 4.488, 4.488, 4.3707, 4.3148]
        )
        expected = np.array(
            [200.3, 200.254, 200.427, 200.735, 200.916, 200.979, 200.984, 200.991, 195.064, 192.356]
        )
        potentialDensity = np.array(
            [-2.7113, -2.6188, -2.5977, -2.6552, -2.6903, -2.7279, -2.7414, -2.7768, 0.6655, 1.7939]
        )
        # fmt: on
        result = dc.convert_oxygen_to_umol_per_kg(oxMlPerL, potentialDensity)
        request.node.return_value = result.tolist()
        assert np.allclose(expected, result, rtol=0, atol=1e-2)

    def test_convert_to_umol_per_l(self, request):
        # fmt: off
        oxMlPerL = np.array(
            [4.4728, 4.4722, 4.4762, 4.4828, 4.4867, 4.4879, 4.488, 4.488, 4.3707, 4.3148]
        )
        expected = np.array(
            [199.757, 199.730, 199.906, 200.202, 200.375, 200.430, 200.433, 200.433, 195.194, 192.701]
        )
        # fmt: on
        result = dc.convert_oxygen_to_umol_per_l(oxMlPerL)
        request.node.return_value = result.tolist()
        assert np.allclose(expected, result, rtol=0, atol=1e-2)

    def test_convert_to_dov_dt(self, request, source_data):
        expected = source_data["sbox1dV/dT"].values
        result = dc.convert_sbe43_oxygen(
            source_data["sbeox0V"].values,
            source_data["tv290C"].values,
            source_data["prdM"].values,
            source_data["sal00"].values,
            tc.oxygen_43_coefs_sn1686,
            False,
            False,
            2,
            0.25,
            "dov/dt",
        )

        request.node.return_value = result.tolist()
        # TODO: Fix this test
        assert np.allclose(expected, result, rtol=0, atol=1e-5)

    def test_convert_to_pct_saturation(self, request, source_data):
        expected = source_data["sbeox0PS"].values
        result = dc.convert_sbe43_oxygen(
            source_data["sbeox0V"].values,
            source_data["tv290C"].values,
            source_data["prdM"].values,
            source_data["sal00"].values,
            tc.oxygen_43_coefs_sn1686,
            False,
            False,
            2,
            0.25,
            "saturation_percent",
        )
        request.node.return_value = result.tolist()
        assert np.allclose(expected, result, rtol=0, atol=1e-2)

    @pytest.mark.parametrize(
        "from_units,to_units",
        [
            ("ml/l", "ml/l"),
            ("ml/l", "mg/l"),
            ("ml/l", "umol/kg"),
            ("ml/l", "umol/l"),
            ("ml/l", "saturation_percent"),
            ("mg/l", "ml/l"),
            ("mg/l", "mg/l"),
            ("mg/l", "umol/kg"),
            ("mg/l", "umol/l"),
            ("mg/l", "saturation_percent"),
            ("umol/kg", "ml/l"),
            ("umol/kg", "mg/l"),
            ("umol/kg", "umol/kg"),
            ("umol/kg", "umol/l"),
            ("umol/kg", "saturation_percent"),
            ("umol/l", "ml/l"),
            ("umol/l", "mg/l"),
            ("umol/l", "umol/kg"),
            ("umol/l", "umol/l"),
            ("umol/l", "saturation_percent"),
            ("saturation_percent", "ml/l"),
            ("saturation_percent", "mg/l"),
            ("saturation_percent", "umol/kg"),
            ("saturation_percent", "umol/l"),
            ("saturation_percent", "saturation_percent"),
        ],
    )
    def test_convert_oxygen_units_all_combinations(self, source_data, from_units, to_units):
        temperature = source_data["tv290C"].values
        pressure = source_data["prdM"].values
        salinity = source_data["sal00"].values

        oxygen_ml_per_l = source_data["sbeox0ML/L"].values
        oxygen_mg_per_l = source_data["sbeox0Mg/L"].values
        oxygen_umol_per_l = source_data["sbeox0Mm/L"].values
        oxygen_umol_per_kg = source_data["sbox0Mm/Kg"].values
        oxygen_saturation_percent = source_data["sbeox0PS"].values

        # Convert to oxygen_in units from source data, otherwise we start from rounded data and lose precision in the conversions.
        oxygen_in = dc.convert_sbe43_oxygen(
            source_data["sbeox0V"].values,
            source_data["tv290C"].values,
            source_data["prdM"].values,
            source_data["sal00"].values,
            tc.oxygen_43_coefs_sn1686,
            False,
            False,
            2,
            0.25,
            from_units,
        )

        # Step 2: convert expected ml/l values to the target units.
        if to_units == "ml/l":
            expected = oxygen_ml_per_l
        elif to_units == "mg/l":
            expected = oxygen_mg_per_l
        elif to_units == "umol/l":
            expected = oxygen_umol_per_l
        elif to_units == "umol/kg":
            expected = oxygen_umol_per_kg
        elif to_units == "saturation_percent":
            expected = oxygen_saturation_percent
        else:
            raise ValueError(f"unsupported to_units in test: {to_units}")

        # Function preserves explicit bad flags from input.
        expected = np.where(oxygen_in == const.FLAG_VALUE, const.FLAG_VALUE, expected)

        result = dc.convert_oxygen_units(
            oxygen_in,
            temperature,
            pressure,
            salinity,
            from_units,
            to_units,
        )

        atol_by_to_units = {
            "ml/l": 1e-3,
            "mg/l": 1e-3,
            "umol/kg": 1e-2,
            "umol/l": 1e-2,
            "saturation_percent": 1e-2,
        }
        assert np.allclose(expected, result, rtol=0, atol=atol_by_to_units[to_units])


class TestConvertChlorophylla:
    def test_convert_eco(self, request):
        # fmt: off
        rawAnalog = np.array(
            [0.0949, 0.0948, 0.0960, 0.0961, 0.0962, 0.0959, 0.1013, 0.1012, 0.1015, 0.1012, 0.1003, 0.0999, 0.0999, 0.0996]
        )
        expected = np.array(
            [0.2691, 0.2683, 0.2798, 0.2813, 0.2821, 0.279, 0.3332, 0.3317, 0.3355, 0.3317, 0.3233, 0.3187, 0.3195, 0.3157]
        )
        # fmt: on
        result = dc.convert_eco(rawAnalog, tc.chlorophyll_a_coefs_sn6130)
        request.node.return_value = result.tolist()
        assert np.allclose(expected, result, rtol=0, atol=1e-2)


class TestConvertpH:
    def test_convert_sbe18_ph(self, request):
        # fmt: off
        rawVolts = np.array(
            [2.9507, 2.9518, 2.9522, 2.9517, 2.9521, 2.9519, 2.9517, 2.9526, 2.9525, 2.9537, 2.953, 2.9527]
        )
        temperatureC = np.array(
            [26.6927, 26.6994, 26.6624, 26.6122, 26.5808, 26.5621, 26.5503, 26.5505, 26.569, 26.5954, 26.6182, 26.6359]
        )
        expected = np.array(
            [8.587, 8.591, 8.593, 8.591, 8.593, 8.592, 8.592, 8.595, 8.594, 8.599, 8.596, 8.595]
        )
        # fmt: on
        result = dc.convert_sbe18_ph(rawVolts, temperatureC, tc.ph_coefs_sn0762)
        request.node.return_value = result.tolist()
        assert np.allclose(expected, result, rtol=0, atol=1e-3)


class TestPARlogarithmic:
    def test_convert_par_logarithmic(self, request):
        # fmt: off
        rawVolts = np.array(
            [1.176241702, 1.175249866, 1.174944686, 1.176394293, 0.089951934, 0.090562294, 0.090638590, 0.091172655, 0.090485999, 0.090943770, 0.090409704, 0.091096360, 0.093537804]
        )
        expected = np.array(
            [0.81605, 0.81394, 0.81330, 0.81637, 0.04817, 0.04824, 0.04825, 0.04832, 0.04823, 0.04829, 0.04822, 0.04831, 0.04862]
        )
        # fmt :on
        result = dc.convert_par_logarithmic(rawVolts, tc.par_coefs_sn0411)
        request.node.return_value = result.tolist()
        assert np.allclose(expected, result, rtol=0, atol=1e-3)


# class TestContourFromTSP:
#     # Note: this class doesn't actually test anything and is only for debug
#     data_path = test_data / 'SBE37SM-RS232_03711722_2015_11_18_subset_derived.asc'
#     @pytest.fixture
#     def data(self):
#         return pd.read_csv(self.data_path)

#     def test_contour_from_t_s_p_pass(self, data, request):
#         temperature = self.data['t090C'].values
#         salinity = self.data['sal00'].values
#         pressure = self.data['prM'].values
#         contour_data = c.contour_from_t_s_p(
#             temperature, salinity, pressure
#         )
#         assert(True)


class TestSeaFETPH:
    internal_ph_counts = np.array([5105334, 5105384, 5105350, 5105505, 5105347])
    external_ph_counts = np.array([5139133, 5139214, 5139120, 5139294, 5139183])
    expected_internal_ph_volts = np.array([-0.978492, -0.978477, -0.978487, -0.978441, -0.978488])
    expected_external_ph_volts = np.array([-0.968419, -0.968395, -0.968423, -0.968371, -0.968404])
    expected_internal_ph = np.array([6.8897, 6.89, 6.8898, 6.8906, 6.8898])
    expected_external_ph = np.array([6.2466, 6.247, 6.2466, 6.2474, 6.2469])
    ph_temperature = np.array([22.751, 22.7514, 22.752, 22.752, 22.752])

    def ph_voltage_to_counts(self, voltage: np.ndarray):
        adc_vref = 2.5
        gain = 1
        adc_23bit = 8388608
        ph_counts = adc_23bit * (voltage * gain / adc_vref + 1)
        return ph_counts

    def test_convert_ph_voltage_counts(self, request):
        internal_ph_volts = dc.convert_ph_voltage_counts(self.internal_ph_counts)
        request.node.return_value = internal_ph_volts.tolist()
        assert np.allclose(internal_ph_volts, self.expected_internal_ph_volts, atol=1e-6)

    def test_convert_internal_seafet_ph(self, request):
        internal_ph = dc.convert_internal_seafet_ph(
            raw_ph=self.internal_ph_counts,
            temperature=self.ph_temperature,
            coefs=tc.ph_seafet_internal_coefs,
        )
        request.node.return_value = internal_ph.tolist()
        assert np.allclose(internal_ph, self.expected_internal_ph, atol=1e-6)

    def test_convert_external_seafet_ph(self, request):
        external_ph = dc.convert_external_seafet_ph(
            raw_ph=self.external_ph_counts,
            temperature=self.ph_temperature,
            salinity=35,
            pressure=0,
            coefs=tc.ph_seafet_external_coefs,
            ph_units="counts",
        )
        request.node.return_value = external_ph.tolist()
        assert np.allclose(external_ph, self.expected_external_ph, atol=1e-6)

    def test_internal_shallow_ph_app_note_99(self, request):
        # from application note 99
        ph_counts = self.ph_voltage_to_counts(-1.010404)
        # coefs = cc.PHSeaFETInternalCoefficients(k0=-1.438788, k2=-1.304895e-3)
        internal_ph = dc.convert_internal_seafet_ph(
            raw_ph=ph_counts,
            temperature=15.8735,
            coefs=tc.internal_ph_coefs,
        )
        assert np.allclose(internal_ph, 7.8310, atol=1e-4)

    def test_external_shallow_seafet_seaphox_ph_app_note_99(self, request):
        # from application note 99
        # Shallow pH is effectively the same as deep pH with pressure
        # set to 0. Since the deep pH example in app note 99 is valid,
        # there must be a typo in the shellow pH example. This test uses
        # the same inputs as the shallow example and fixes the result
        external_ph = dc.convert_external_seafet_ph(
            raw_ph=-0.965858,
            temperature=15.8735,
            salinity=36.817,
            pressure=0,
            coefs=tc.external_shallow_ph_coefs,
            ph_units="volts",
        )
        assert np.allclose(external_ph, 7.9224, atol=1e-4)

    def test_external_deep_seaphox_float_ph_app_note_99(self, request):
        # from application note 99
        external_ph = dc.convert_external_seafet_ph(
            raw_ph=-0.885081,
            temperature=23.4169,
            salinity=34.812,
            pressure=100,
            coefs=tc.external_ph_coefs,
            ph_units="volts",
        )
        assert np.allclose(external_ph, 7.9394, atol=1e-4)

    def test_convert_external_seafet_ph_legacy(self, request):
        # from application note 99
        external_ph = dc.convert_external_seafet_ph(
            raw_ph=self.ph_voltage_to_counts(-0.885081),
            temperature=23.4169,
            salinity=34.812,
            pressure=100,
            coefs=tc.external_ph_coefs,
            ph_units="counts",
            formula_version="legacy",
        )
        assert np.allclose(external_ph, 7.939415, atol=1e-6)

    def test_convert_external_seafet_ph_legacy_volts(self, request):
        # from application note 99
        external_ph = dc.convert_external_seafet_ph(
            raw_ph=-0.885081,
            temperature=23.4169,
            salinity=34.812,
            pressure=100,
            coefs=tc.external_ph_coefs,
            ph_units="volts",
            formula_version="legacy",
        )
        assert np.allclose(external_ph, 7.939415, atol=1e-6)

    def test_convert_external_seafet_ph_v1p3(self, request):
        # from application note 99
        external_ph = dc.convert_external_seafet_ph(
            raw_ph=-0.885081,
            temperature=23.4169,
            salinity=34.812,
            pressure=100,
            coefs=tc.external_ph_coefs,
            ph_units="volts",
        )
        assert np.allclose(external_ph, 7.939406, atol=1e-6)

    def test_convert_external_seafet_ph_v1p3_shallow(self, request):
        # from application note 99
        external_ph = dc.convert_external_seafet_ph(
            raw_ph=-0.900081,
            temperature=23.4169,
            salinity=34.812,
            pressure=0,
            coefs=tc.external_ph_coefs,
            ph_units="volts",
        )
        assert np.allclose(external_ph, 7.6653, atol=1e-4)

    def test_convert_external_seafet_ph_argo_example(self, request):
        # from Processing BGC-Argo pH data at the DAC level
        # https://archimer.ifremer.fr/doc/00460/57195/
        expected_ph = [7.9870, 7.9788, 7.9105, 7.7846, 7.7326, 7.7307]
        external_ph = dc.convert_external_seafet_ph(
            raw_ph=np.array([-0.82966, -0.83049, -0.83921, -0.85150, -0.85968, -0.86191]),
            temperature=np.array([29.545, 29.607, 25.741, 21.094, 15.725, 13.407]),
            salinity=np.array([34.637, 34.767, 35.104, 35.080, 34.886, 34.890]),
            pressure=np.array([2.2, 42, 82, 122, 162, 202]),
            coefs=tc.external_ph_k2_poly_coefs,
            ph_units="volts",
            formula_version="1.3",
        )
        assert np.allclose(external_ph, expected_ph, atol=1e-4)


class TestInternalSeaFETTemperature:
    def test_convert_internal_seafet_temperature(self, request):
        temperature_counts = np.array([25616, 25600])
        expected_temperature = np.array([21.8335, 21.7906])
        temperature = dc.convert_internal_seafet_temperature(temperature_counts=temperature_counts)
        request.node.return_value = temperature.tolist()
        assert np.allclose(expected_temperature, temperature, atol=1e-6)


class TestSeaFETRelativeHumidity:
    def test_convert_seafet_relative_humidity(self, request):
        humidity_counts = np.array([24096, 24160])
        temperature_counts = np.array([25616, 25600])
        expected_humidity = np.array([39.4845, 39.6001])

        temperature = dc.convert_internal_seafet_temperature(temperature_counts=temperature_counts)
        humidity = dc.convert_seafet_relative_humidity(
            humidity_counts=humidity_counts, temperature=temperature
        )

        request.node.return_value = humidity.tolist()
        assert np.allclose(expected_humidity, humidity, atol=1e-6)


class TestConvertSBE63Oxygen:
    raw_oxygen = np.array([31.06, 31.66, 32.59, 33.92, 34.82, 35.44, 35.44, 35.44])
    pressure = np.array([0, 0, 0, 0, 0, 0, 1000, 100])
    temperature = np.array([30, 26, 20, 12, 6, 2, 2, 2])
    salinity = np.array([0, 0, 0, 0, 0, 0, 0, 35])  #  salinity is 0 PSU during calibration
    expected_oxygen = np.array([0.706, 0.74, 0.799, 0.892, 1.005, 1.095, 1.1398, 0.8647])

    # Some tests need to be run on complete datasets
    cnv_path = test_data / "SBE37SMP-ODO-RS232_03711459_2023_08_11-seasoft-convert-o2-full.cnv"

    @pytest.fixture
    def source_data(self):
        return id.read_cnv_file(self.cnv_path, "seasoft")

    def test_convert_sbe63_oxygen(self, request):
        oxygen = dc.convert_sbe63_oxygen(
            self.raw_oxygen,
            self.temperature,
            self.pressure,
            self.salinity,
            tc.oxygen_63_coefs_sn2568,
            tc.thermistor_63_coefs_sn2568,
            "C",
        )
        # TODO: fix tolerance in TKIT-110
        request.node.return_value = oxygen.tolist()
        assert np.allclose(self.expected_oxygen, oxygen, atol=1e-3)

    def test_thermistor_temperature(self, request):
        raw_temperature = np.array([1.12015, 1.12015, 1.12016, 1.12016, 0.99562, 0.82934, 0.64528])

        expected = np.array([2.0002, 2.0002, 1.9999, 1.9999, 5.9998, 12.0, 19.9998])
        thermistor_temperature = dc.convert_sbe63_thermistor(
            raw_temperature, tc.thermistor_63_coefs_sn2568
        )

        request.node.return_value = thermistor_temperature.tolist()
        assert np.allclose(expected, thermistor_temperature, atol=1e-6)

    def test_convert_sbe63_oxygen_therm_volts(self, request):
        raw_oxygen = np.array([31.06, 31.66, 32.59, 33.92, 34.82, 35.44])
        pressure = np.array([0, 0, 0, 0, 0, 0])
        raw_temperature = np.array([0.6, 0.5, 0.4, 0.35, 0.3, 0.25])
        salinity = np.array([0, 0, 0, 0, 0, 0])
        expected = np.array([0.93, 0.688, 0.459, 0.304, 0.206, 0.137])

        result = dc.convert_sbe63_oxygen(
            raw_oxygen,
            raw_temperature,
            pressure,
            salinity,
            tc.oxygen_63_coefs_sn2568,
            tc.thermistor_63_coefs_sn2568,
        )
        request.node.return_value = result.tolist()
        assert np.allclose(expected, result, rtol=0, atol=1e-3)

    def test_convert_sbe63_oxygen_from_hex(self, request):
        raw_oxygen = np.array([16.774, 16.775, 16.779, 16.778, 16.774, 16.779])
        pressure = np.array([-0.057, -0.062, -0.057, -0.056, -0.068, -0.056])
        raw_temperature = np.array([0.581763, 0.581758, 0.581753, 0.581737, 0.581725, 0.581717])
        salinity = np.array([0.0115, 0.0115, 0.0115, 0.0115, 0.0115, 0.0115])
        expected = np.array([5.872, 5.872, 5.868, 5.869, 5.872, 5.868])

        result = dc.convert_sbe63_oxygen(
            raw_oxygen,
            raw_temperature,
            pressure,
            salinity,
            tc.oxygen_63_coefs_sn11459,
            tc.thermistor_63_coefs_sn11459,
        )
        request.node.return_value = result.tolist()
        assert np.allclose(expected, result, rtol=0, atol=1e-3)

    def test_convert_sbe63_oxygen_ml_per_l(self, request, source_data):
        expected = source_data["sbeopoxML/L"].values
        result = dc.convert_sbe63_oxygen(
            source_data["sbeoxpd"].values,
            source_data["sbeoxtv"].values,
            source_data["prdM"].values,
            source_data["sal00"].values,
            tc.oxygen_63_coefs_sn11459,
            tc.thermistor_63_coefs_sn11459,
            "volts",
            "ml/l",
        )
        request.node.return_value = result.tolist()
        assert np.allclose(expected, result, rtol=0, atol=1e-3)

    def test_convert_sbe63_oxygen_mg_per_l(self, request, source_data):
        expected = source_data["sbeopoxMg/L"].values
        result = dc.convert_sbe63_oxygen(
            source_data["sbeoxpd"].values,
            source_data["sbeoxtv"].values,
            source_data["prdM"].values,
            source_data["sal00"].values,
            tc.oxygen_63_coefs_sn11459,
            tc.thermistor_63_coefs_sn11459,
            "volts",
            "mg/l",
        )
        request.node.return_value = result.tolist()
        assert np.allclose(expected, result, rtol=0, atol=1e-3)

    def test_convert_sbe63_oxygen_saturation_percent(self, request, source_data):
        expected = source_data["sbeopoxPS"].values
        result = dc.convert_sbe63_oxygen(
            source_data["sbeoxpd"].values,
            source_data["sbeoxtv"].values,
            source_data["prdM"].values,
            source_data["sal00"].values,
            tc.oxygen_63_coefs_sn11459,
            tc.thermistor_63_coefs_sn11459,
            "volts",
            "saturation_percent",
            source_data["tv290C"].values,
        )
        request.node.return_value = result.tolist()
        print(f"Expected: {expected}")
        assert np.allclose(expected, result, rtol=0, atol=1e-2)

    def test_convert_sbe63_oxygen_umol_per_kg(self, request, source_data):
        expected = source_data["sbeopoxMm/Kg"].values
        result = dc.convert_sbe63_oxygen(
            source_data["sbeoxpd"].values,
            source_data["sbeoxtv"].values,
            source_data["prdM"].values,
            source_data["sal00"].values,
            tc.oxygen_63_coefs_sn11459,
            tc.thermistor_63_coefs_sn11459,
            "volts",
            "umol/kg",
            source_data["tv290C"].values,
        )
        request.node.return_value = result.tolist()
        assert np.allclose(expected, result, rtol=0, atol=1e-2)

    def test_convert_sbe63_oxygen_umol_per_l(self, request, source_data):
        expected = source_data["sbeopoxMm/L"].values
        result = dc.convert_sbe63_oxygen(
            source_data["sbeoxpd"].values,
            source_data["sbeoxtv"].values,
            source_data["prdM"].values,
            source_data["sal00"].values,
            tc.oxygen_63_coefs_sn11459,
            tc.thermistor_63_coefs_sn11459,
            "volts",
            "umol/l",
        )
        request.node.return_value = result.tolist()
        assert np.allclose(expected, result, rtol=0, atol=1e-2)

    def test_convert_sbe63_oxygen_raw_pd(self, request, source_data):
        expected = source_data["sbeoxpd"].values
        result = dc.convert_sbe63_oxygen(
            source_data["sbeoxpd"].values,
            source_data["sbeoxtv"].values,
            source_data["prdM"].values,
            source_data["sal00"].values,
            tc.oxygen_63_coefs_sn11459,
            tc.thermistor_63_coefs_sn11459,
            "volts",
            "raw_phase_usec",
        )
        request.node.return_value = result.tolist()
        assert np.allclose(expected, result, rtol=0, atol=1e-2)

    def test_convert_sbe63_oxygen_raw_pd_v(self, request, source_data):
        expected = source_data["sbeoxpdv"].values
        result = dc.convert_sbe63_oxygen(
            source_data["sbeoxpd"].values,
            source_data["sbeoxtv"].values,
            source_data["prdM"].values,
            source_data["sal00"].values,
            tc.oxygen_63_coefs_sn11459,
            tc.thermistor_63_coefs_sn11459,
            "volts",
            "raw_phase_v",
        )
        request.node.return_value = result.tolist()
        assert np.allclose(expected, result, rtol=0, atol=1e-3)

    def test_convert_sbe63_oxygen_ox_temp_c(self, request, source_data):
        expected = source_data["sbeoxTC"].values
        result = dc.convert_sbe63_oxygen(
            source_data["sbeoxpd"].values,
            source_data["sbeoxtv"].values,
            source_data["prdM"].values,
            source_data["sal00"].values,
            tc.oxygen_63_coefs_sn11459,
            tc.thermistor_63_coefs_sn11459,
            "volts",
            "ox_temperature_c",
        )
        request.node.return_value = result.tolist()
        assert np.allclose(expected, result, rtol=0, atol=1e-3)

    def test_convert_sbe63_oxygen_ox_temp_f(self, request, source_data):
        expected = source_data["sbeoxTF"].values
        result = dc.convert_sbe63_oxygen(
            source_data["sbeoxpd"].values,
            source_data["sbeoxtv"].values,
            source_data["prdM"].values,
            source_data["sal00"].values,
            tc.oxygen_63_coefs_sn11459,
            tc.thermistor_63_coefs_sn11459,
            "volts",
            "ox_temperature_f",
        )
        request.node.return_value = result.tolist()
        assert np.allclose(expected, result, rtol=0, atol=1e-3)

    def test_oxygen_saturation_percent_gg(self, request, source_data):
        result = dc.derive_oxygen_saturation_gg(
            source_data["tv290C"].values, source_data["sal00"].values
        )
        request.node.return_value = result.tolist()
        assert np.allclose(result, source_data["oxsolML/L"].values, rtol=0, atol=1e-4)


class TestSPAR:
    def test_surface_par(self, request):
        spar_raw = np.array([0.03296703296703297])
        expected_spar_biospherical = np.array([51.41538])  # from SBE data processing
        expected_spar_linear = np.array([-6636])  # from SBE data processing
        expected_spar_logarithmic = np.array([1127.0])  # from SBE data processing

        spar_biospherical = dc.convert_spar_biospherical(spar_raw, tc.spar_coefs)
        spar_linear = dc.convert_spar_linear(spar_raw, tc.spar_coefs)
        spar_logarithmic = dc.convert_spar_logarithmic(spar_raw, tc.spar_coefs)

        assert np.allclose(expected_spar_biospherical, spar_biospherical, atol=1e-5)
        assert np.allclose(expected_spar_linear, spar_linear, atol=0)
        assert np.allclose(expected_spar_logarithmic, spar_logarithmic, atol=1e-1)


class TestAltimeter:
    def test_altimeter(self, request):
        altimeter_raw = np.array([3.1893, 1.1807, 2.3797])
        expected = np.array([63.79, 23.61, 47.59])  # meters, from SBE data processing

        height = dc.convert_altimeter(altimeter_raw, tc.altimeter_coefs)

        assert np.allclose(expected, height, atol=1e-2)


class TestBuoyancy:
    # fmt: off
    # Testing data comes from a CalCOFI cruise
    # SBE911plus\Fathom-Testing\RL2301001seasoft-convert-bin-5-buoy2.cnv
    temperature = np.asarray([13.4288, 10.3983, 9.2891, 8.3246, 7.6800, 7.1504, 6.7090, 6.1575, 5.8453, 5.6158])
    # bin_average will return floats so make sure we're replicating that here
    pressure = np.asarray([50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0])
    salinity = np.array([33.1678, 33.6430, 33.9238, 33.9917, 34.0375, 34.0499, 34.0697, 34.1113, 34.1726, 34.2101])

    expected_buoyancy_freq_squared = np.asarray([-9.990e-29, 1.3072e-04, 5.9338e-05, 3.3027e-05, 2.1500e-05, 1.6153e-05, 1.8199e-05, 1.9350e-05, 1.4610e-05, -9.990e-29])
    expected_buoyancy_freq = np.asarray([-9.990e-29, 6.55, 4.41, 3.29, 2.66, 2.30, 2.44, 2.52, 2.19, -9.990e-29])
    expected_stability = np.asarray([-9.990e-29, 1.3343e-05, 6.0570e-06, 3.3713e-06, 2.1946e-06, 1.6489e-06, 1.8577e-06, 1.9752e-06, 1.4913e-06, -9.990e-29])
    expected_scaled_stability = np.asarray([-9.990e-29, 1334.3, 605.7, 337.1, 219.5, 164.9, 185.8, 197.5, 149.1, -9.990e-29])
    # fmt: on

    def test_buoyancy(self):
        (buoyancy_freq_squared, buoyancy_freq, stability, scaled_stability) = dc.buoyancy(
            self.temperature,
            self.salinity,
            self.pressure,
            np.asarray([34.034167]),  # converted from metadata 34.02.03 N in H,M,S
            np.asarray([121.060556]),  # converted from metadata 121 03.38 W in H, M, S
            150,  # window size
            True,
        )

        # Comparing EOS-80 to TEOS-10 buoyancy calculations.
        # We do not expect them to agree better than +/-1.5% due to differences in the algorithms
        rel_tol = 0.02  # 2%
        assert buoyancy_freq_squared == pytest.approx(
            self.expected_buoyancy_freq_squared, rel=rel_tol
        )
        assert buoyancy_freq == pytest.approx(self.expected_buoyancy_freq, rel=rel_tol)
        assert stability == pytest.approx(self.expected_stability, rel=rel_tol)
        assert scaled_stability == pytest.approx(self.expected_scaled_stability, rel=rel_tol)

        # fmt: off
        # adding exact result comparisons to detect changes that still pass the tolerance tests
        expected_buoyancy_freq_squared = np.array([-9.99e-29, 0.00012830774423210731, 5.9061904158766353e-05, 3.285253714504237e-05, 2.1427260827130063e-05, 1.610786846795634e-05, 1.8197855102567525e-05, 1.928895593187585e-05, 1.4586979090250307e-05, -9.99e-29])
        expected_buoyancy_freq = np.array([-9.99e-29, 6.490065311850448, 4.403280527245377, 3.28402980427617, 2.6521981054713932, 2.2995437132561976, 2.4441774544272, 2.5163844503226995, 2.188292201351598, -9.99e-29])
        expected_scaled_stability = np.array([-9.99e-29, 1309.6979919526298, 602.8661047814034, 335.33390726433163, 218.7108575931272, 164.41328104617097, 185.74372883265428, 196.8782837049103, 148.88453614380862, -9.99e-29])
        expected_stability = np.array([-9.99e-29, 1.3096979919526299e-05, 6.028661047814034e-06, 3.353339072643316e-06, 2.187108575931272e-06, 1.6441328104617097e-06, 1.8574372883265428e-06, 1.968782837049103e-06, 1.4888453614380862e-06, -9.99e-29])
        # fmt: on
        assert np.all(buoyancy_freq_squared == expected_buoyancy_freq_squared)
        assert np.all(buoyancy_freq == expected_buoyancy_freq)
        assert np.all(stability == expected_stability)
        assert np.all(scaled_stability == expected_scaled_stability)

    def test_buoyancy_eos80(self):
        (buoyancy_freq_squared, buoyancy_freq, stability, scaled_stability) = dc.buoyancy(
            self.temperature,
            self.salinity,
            self.pressure,
            np.asarray([34.034167]),  # converted from metadata 34.02.03 N in H,M,S
            np.asarray([121.060556]),  # converted from metadata 121 03.38 W in H, M, S
            150,  # window size
            False,
        )

        # Comparing SBE Data Processing C++ to local Python results using the same EOS-80 calculations.
        # We expect very very close agreement: << 1% differnce
        rel_tol = 0.002  # 0.2%
        assert buoyancy_freq_squared == pytest.approx(
            self.expected_buoyancy_freq_squared, rel=rel_tol
        )
        assert buoyancy_freq == pytest.approx(self.expected_buoyancy_freq, rel=rel_tol)
        assert stability == pytest.approx(self.expected_stability, rel=rel_tol)
        assert scaled_stability == pytest.approx(self.expected_scaled_stability, rel=rel_tol)

        # fmt: off
        # adding exact result comparisons to detect changes that still pass the tolerance tests
        expected_buoyancy_freq_squared = np.array([-9.99e-29, 0.00013069443519014734, 5.932961410753015e-05, 3.3021309348275255e-05, 2.1495844621369152e-05, 1.6150523253891014e-05, 1.8197348677595506e-05, 1.9348996910895833e-05, 1.4610475131109874e-05, -9.99e-29])
        expected_buoyancy_freq = np.array([-9.99e-29, 6.550149019323256, 4.4132486213213316, 3.2924544645935523, 2.6564392562582047, 2.302586378268702, 2.4441434448941073, 2.520297798385073, 2.1900538928778936, -9.99e-29])
        expected_scaled_stability = np.array([-9.99e-29, 1334.060078386804, 605.5987165439468, 337.05660655229923, 219.41090136335237, 164.84865915481353, 185.73855980245008, 197.49110924835415, 149.12435256658875, -9.99e-29])
        expected_stability = np.array([-9.99e-29, 1.334060078386804e-05, 6.055987165439468e-06, 3.370566065522992e-06, 2.1941090136335238e-06, 1.6484865915481353e-06, 1.857385598024501e-06, 1.9749110924835414e-06, 1.4912435256658876e-06, -9.99e-29])
        # fmt: on
        assert np.all(buoyancy_freq_squared == expected_buoyancy_freq_squared)
        assert np.all(buoyancy_freq == expected_buoyancy_freq)
        assert np.all(stability == expected_stability)
        assert np.all(scaled_stability == expected_scaled_stability)


class TestDeriveDescentRateAcceleration:
    cnv_path = test_data / "SBE19plus_01906398_2019_07_15_0033-seasoft-convert-speeds.cnv"

    @pytest.fixture
    def source_data(self):
        return id.read_cnv_file(self.cnv_path, "seasoft")

    def test_derive_descent_rate_meters(self, source_data):
        descent_rate_m = dc.derive_descent_rate(source_data["depSM"].values, 2, 0.25)
        assert np.allclose(descent_rate_m, source_data["dz/dtM"].values, rtol=0, atol=1e-2)

    def test_derive_descent_rate_feet(self, source_data):
        descent_rate_f = dc.derive_descent_rate(source_data["depSF"].values, 2, 0.25)
        # TODO: SBE data processing is imprecise and returns a value that is slightly different than the expected value. The atol is set to 1e-1 to account for this.
        expected_dzdtF = (
            source_data["dz/dtF"].values * 3.28084
        )  # TODO: for some reason the dz/dtF returns meters/s too in SBE Data Proc
        assert np.allclose(descent_rate_f, expected_dzdtF, rtol=0, atol=1e-1)

    def test_derive_acc_meters(self, source_data):
        acc_m = dc.derive_acceleration(source_data["depSM"].values, 2, 0.25)
        assert np.allclose(acc_m, source_data["accM"].values, rtol=0, atol=1e-2)

    def test_derive_acc_feet(self, source_data):
        acc_f = dc.derive_acceleration(source_data["depSF"].values, 2, 0.25)
        # TODO: SBE data processing is imprecise and returns a value that is slightly different than the expected value. The atol is set to 1e-1 to account for this.
        assert np.allclose(acc_f, source_data["accF"].values, rtol=0, atol=1e-1)


class TestDeriveOxygenSaturation:
    cnv_path = test_data / "SBE19plus_01906398_2019_07_15_0033-seasoft-convert-Ox-Sat.cnv"

    @pytest.fixture
    def source_data(self):
        return id.read_cnv_file(self.cnv_path, "seasoft")

    def test_derive_oxygen_saturation_gg(self, source_data):
        ox_sat_gg = dc.derive_oxygen_saturation_gg(
            source_data["tv290C"].values, source_data["sal00"].values
        )
        assert np.allclose(ox_sat_gg, source_data["oxsolML/L"].values, rtol=0, atol=1e-3)

    def test_derive_oxygen_saturation_w(self, source_data):
        ox_sat_gg = dc.derive_oxygen_saturation_w(
            source_data["tv290C"].values, source_data["sal00"].values
        )
        assert np.allclose(ox_sat_gg, source_data["oxsatML/L"].values, rtol=0, atol=1e-3)


class TestCstar:
    def test_cstar_attenuation(self, request):
        raw = np.array([3.7186, 4.3632, 4.2264])
        expected = np.array([0.9273, 0.2869, 0.4145])  # from SBE data processing

        result = dc.convert_cstar_attenuation(raw, tc.CSTAR_coefs)

        request.node.return_value = result.tolist()
        assert np.allclose(expected, result, rtol=0, atol=1e-4)

    def test_cstar_transmittance(self, request):
        raw = np.array([3.7186, 4.3632, 4.2264])
        expected = np.array([79.3092, 93.0779, 90.1573])  # from SBE data processing

        result = dc.convert_cstar_transmittance(raw, tc.CSTAR_coefs)

        request.node.return_value = result.tolist()
        assert np.allclose(expected, result, rtol=0, atol=1e-3)


class TestDeriveThermostericAnomaly:
    cnv_path = test_data / "SBE19plus_derive_testing.cnv"

    @pytest.fixture
    def source_data(self):
        return id.read_cnv_file(self.cnv_path, "seasoft")

    def test_derive_tsa(self, source_data):
        tsa = eos80.derive_thermosteric_anomaly(
            source_data["sal00"].values, source_data["tv290C"].values
        )
        assert np.allclose(tsa, source_data["tsa"].values, rtol=0, atol=1e-2)


class TestDeriveSpecificConductance:
    cnv_path = test_data / "SBE19plus_derive_testing.cnv"

    @pytest.fixture
    def source_data(self):
        return id.read_cnv_file(self.cnv_path, "seasoft")

    def test_derive_sc(self, source_data):
        # TODO: no explanation for why this is so off?
        tsa = eos80.derive_specific_conductance(
            source_data["tv290C"].values, source_data["c0S/m"].values, to_units="uS/cm"
        )

        differences = tsa - source_data["specc"].values
        print("\nDifferences between tsa and specc:")
        print(differences)
        print("\nMax difference:", np.max(np.abs(differences)))

        # These use large numbers, use rtol=1e-2 to allow for small relative differences
        assert np.allclose(tsa, source_data["specc"].values, rtol=1e-2, atol=0)


class TestDerivePotentialTemperatureAnomaly:
    cnv_path = test_data / "SBE19plus_derive_testing.cnv"

    @pytest.fixture
    def source_data(self):
        return id.read_cnv_file(self.cnv_path, "seasoft")

    def test_derive_pta(self, source_data):
        pta = eos80.derive_potential_temperature_anomaly(
            source_data["sal00"].values,
            source_data["tv290C"].values,
            source_data["prdM"].values,
            a0=1,
            a1=2,
        )

        assert np.allclose(pta, source_data["pta090C"].values, rtol=0, atol=1e-3)


class TestDeriveGeopotentialAnomaly:
    cnv_path = test_data / "SBE19plus_derive_testing.cnv"

    @pytest.fixture
    def source_data(self):
        return id.read_cnv_file(self.cnv_path, "seasoft")

    def test_derive_gpa(self, source_data):
        # sva = sw.svan(source_data["sal00"].to_numpy(), source_data["tv290C"].to_numpy(), source_data["prdM"].to_numpy()) * 1e8
        gpa = eos80.derive_gpa(source_data["sva"], source_data["prdM"].values)
        differences = gpa - source_data["gpa"].values
        print("\nDifferences between gpa and source:")
        print(differences)
        print("\nMax difference:", np.max(np.abs(differences)))
        assert np.allclose(gpa, source_data["gpa"].values, rtol=0, atol=1e-3)
