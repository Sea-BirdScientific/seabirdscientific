"""Data conversion unit tests.

"""

# Native imports

# Third-party imports
import gsw
import numpy as np
import pandas as pd
import pytest

# Sea-Bird imports

# Internal imports
# import seabirdscientific.cal_coefficients as cc
import seabirdscientific.conversion as dc
import seabirdscientific.instrument_data as id
import seabirdscientific.contour as c

import example_coefficients as ec

# DATASETS = [
#     # cal
#     # cnv
#     # hex
#     # cond_label
#     # use_MV_R
#     (
#         cc.SN6130,
#         id.cnv_to_instrument_data("./tests/resources/test-data/19plus_V2_CTD-processing_example.cnv"),
#         id.read_hex_file("./documentation/example_data/19plus_V2_CTD-processing_example.hex",
#             id.InstrumentType.SBE19Plus,
#             [
#                 id.Sensors.Temperature,
#                 id.Sensors.Conductivity,
#                 id.Sensors.Pressure,
#                 id.Sensors.ExtVolt0,
#                 id.Sensors.ExtVolt1,
#                 id.Sensors.ExtVolt2,
#                 id.Sensors.ExtVolt4,
#             ],
#         ),
#         "c0S/m",
#         True,
#     ),
#     # removing until we get a smaller version of this source (NSI-3079)
#     # (
#     #     cc.SN03706385,
#     #     id.cnv_to_instrument_data("./tests/resources/test-data/SBE37SM.cnv"),
#     #     id.read_hex_file(
#     #         "./tests/resources/test-data/SBE37SM.hex",
#     #         id.InstrumentType.SBE37SM,
#     #         [
#     #             id.Sensors.Temperature,
#     #             id.Sensors.Conductivity,
#     #             id.Sensors.Pressure,
#     #         ],
#     #     ),
#     #     "cond0S/m",
#     #     False,
#     # ),
#     (
#         cc.SN03716125,
#         id.cnv_to_instrument_data("./tests/resources/test-data/SBE37SM-RS232_03716125_2017_11_16.cnv"
#         ),
#         id.read_hex_file("./tests/resources/test-data/SBE37SM-RS232_03716125_2017_11_16.hex",
#             id.InstrumentType.SBE37SM,
#             [
#                 id.Sensors.Temperature,
#                 id.Sensors.Conductivity,
#                 id.Sensors.Pressure,
#             ],
#         ),
#         "cond0S/m",
#         False,
#     ),
# ]


class TestConvertTemperature:
    # test temperature raw values
    test_temp_vals = np.array([322798, 322808, 322827, 322838])

    def test_convert_temperature_90C(self):
        expected = [20.4459, 20.4451, 20.4436, 20.4427]
        result = dc.convert_temperature(
            self.test_temp_vals,
            ec.temperature_coefs_sn6130,
            "ITS90",
            "C",
            True,
        )
        assert np.allclose(expected, result, rtol=0, atol=1e-4)

    def test_convert_temperature_90F(self):
        expected = [68.8027, 68.8012, 68.7984, 68.7968]
        result = dc.convert_temperature(
            self.test_temp_vals,
            ec.temperature_coefs_sn6130,
            "ITS90",
            "F",
            True,
        )
        assert np.allclose(expected, result, rtol=0, atol=1e-4)

    def test_convert_temperature_68C(self):
        expected = [20.4508, 20.4500, 20.4485, 20.4476]
        result = dc.convert_temperature(
            self.test_temp_vals,
            ec.temperature_coefs_sn6130,
            "IPTS68",
            "C",
            True,
        )
        assert np.allclose(expected, result, rtol=0, atol=1e-4)

    def test_convert_temperature_68F(self):
        expected = [68.8115, 68.8100, 68.8073, 68.8057]
        result = dc.convert_temperature(
            self.test_temp_vals,
            ec.temperature_coefs_sn6130,
            "IPTS68",
            "F",
            True,
        )
        assert np.allclose(expected, result, rtol=0, atol=1e-4)


class TestConvertPressure:
    # test pressure raw values
    test_pressure_vals = np.array([533539, 533538, 533540, 533537])

    # test pressure temperature compensation raw values
    test_compensation_vals = np.array([20625, 20626, 20626, 20626]) / 13107

    def test_convert_pressure_psia(self):
        expected = [-0.153, -0.154, -0.151, -0.156]
        result = dc.convert_pressure(
            self.test_pressure_vals,
            self.test_compensation_vals,
            ec.pressure_coefs_sn6130,
            "psia",
        )
        assert np.allclose(expected, result, rtol=0, atol=1e-3)

    def test_convert_pressure_dbar(self):
        expected = [-0.105, -0.106, -0.104, -0.107]
        result = dc.convert_pressure(
            self.test_pressure_vals,
            self.test_compensation_vals,
            ec.pressure_coefs_sn6130,
            "dbar",
        )
        assert np.allclose(expected, result, rtol=0, atol=1e-3)


class TestConductivity19plus:
    cnv_path = "./documentation/example_data/19plus_V2_CTD-processing_example.hex"
    expected_data = id.cnv_to_instrument_data(cnv_path)

    hex_path = "./documentation/example_data/19plus_V2_CTD-processing_example.hex"
    raw = id.read_hex_file(
        hex_path,
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
        ec.temperature_coefs_sn6130,
        "IPTS68",
        "C",
        True,
    )

    # test_press_data = np.array([-0.105, -0.106, -0.104, -0.107, -0.105, -0.104])
    pressure = dc.convert_pressure(
        raw["pressure"][0:6].values,
        raw["temperature compensation"][0:6].values,
        ec.pressure_coefs_sn6130,
        "psia",
    )

    def test_convert_conductivity(self):
        expected = [0.008453, 0.008420, 0.008423, 0.008402, 0.008380, 0.008344]
        result = dc.convert_conductivity(
            self.raw["conductivity"][0:6].values,
            self.temperature,
            self.pressure,
            ec.conductivity_coefs_sn6130,
        )
        assert np.allclose(expected, result, rtol=0, atol=1e-6)


class TestConductivity37SM:
    cnv_path = "./tests/resources/test-data/SBE37SM-RS232_03716125_2017_11_16.cnv"
    expected_data = id.cnv_to_instrument_data(cnv_path)

    hex_path = "./tests/resources/test-data/SBE37SM-RS232_03716125_2017_11_16.hex"
    raw = id.read_hex_file(
        hex_path,
        id.InstrumentType.SBE37SM,
        [
            id.Sensors.Temperature,
            id.Sensors.Conductivity,
            id.Sensors.Pressure,
        ],
    )

    temperature = dc.convert_temperature(
        raw["temperature"][0:6],
        ec.temperature_coefs_sn6130,
        "IPTS68",
        "C",
        True,
    )

    # test_press_data = np.array([-0.105, -0.106, -0.104, -0.107, -0.105, -0.104])
    pressure = dc.convert_pressure(
        raw["pressure"][0:6],
        raw["temperature compensation"][0:6],
        ec.pressure_coefs_sn16125,
        "psia",
    )

    def test_convert_conductivity(self):
        expected = [2.711842, 2.715786, 2.715857, 2.715846, 2.715846, 2.715857]
        result = dc.convert_conductivity(
            self.raw["conductivity"][0:6],
            self.temperature,
            self.pressure,
            ec.conductivity_coefs_sn16125,
        )
        assert np.allclose(expected, result, rtol=0, atol=1e-4)


class TestDeriveDensity:
    data_path = "./tests/resources/test-data/SBE37SM-derived.asc"
    data = pd.read_csv(data_path)

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
    def test_derive_potential_density_from_t_s_p_pass(self, reference_pressure, expected_column):
        temperature = self.data["t090C"].values
        salinity = self.data["sal00"].values
        pressure = self.data["prM"].values
        expected = self.data[expected_column].values
        result = dc.potential_density_from_t_s_p(
            temperature,
            salinity,
            pressure,
            reference_pressure=reference_pressure,
        )
        # TODO: improve passing condition for all instances of allclose
        assert np.allclose(result, expected, atol=0.1)

    def test_derive_density_from_t_s_p_pass(self):
        temperature = self.data["t090C"].values
        salinity = self.data["sal00"].values
        pressure = self.data["prM"].values
        expected = self.data["gsw_densityA0"].values
        result = dc.density_from_t_s_p(temperature, salinity, pressure)

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
    def test_derive_potential_density_from_t_c_p_pass(self, reference_pressure, expected_column):
        temperature = self.data["t090C"].values
        conductivity = self.data["c0S/m"].values * 10.0
        pressure = self.data["prM"].values
        expected = self.data[expected_column].values
        result = dc.potential_density_from_t_c_p(
            temperature,
            conductivity,
            pressure,
            reference_pressure=reference_pressure,
        )
        assert np.allclose(result, expected, atol=0.1)

    def test_derive_density_from_t_c_p_pass(self):
        temperature = self.data["t090C"].values
        conductivity = self.data["c0S/m"].values * 10.0
        pressure = self.data["prM"].values
        expected = self.data["gsw_densityA0"].values
        result = dc.density_from_t_c_p(temperature, conductivity, pressure)

        assert np.allclose(result, expected, atol=0.1)


class TestDepthFromPressure:
    data_path = "./tests/resources/test-data/SBE37SM.asc"
    data = pd.read_csv(data_path).loc[6:10]

    @pytest.mark.parametrize(
        "pressure, pressure_units, depth, depth_units",
        [
            ("prdM", "dbar", "depSM", "m"),
            ("prdM", "dbar", "depSF", "ft"),
            ("prdE", "psi", "depSM", "m"),
            ("prdE", "psi", "depSF", "ft"),
        ],
    )
    def test_depth_from_pressure_pass(self, pressure, pressure_units, depth, depth_units):
        expected_depth = self.data[depth].values
        pressure = self.data[pressure].values
        result_depth = dc.depth_from_pressure(pressure, 0, depth_units, pressure_units)
        assert np.allclose(expected_depth, result_depth, atol=0.002)


# class TestSalinityFromTCP:
#     @pytest.mark.parametrize("cal, cnv, hex, cond_label, use_MV_R", DATASETS)
#     def test_salinity_from_tcp(self, cal, cnv, hex, cond_label, use_MV_R):
#         expected_salinity = cnv.measurements["sal00"].values[800:810]
#         temperature = cnv.measurements["tv290C"].values[800:810]
#         pressure = cnv.measurements["prdM"].values[800:810]
#         conductivity = cnv.measurements[cond_label].values[800:810] * 10
#         result_salinity = gsw.SP_from_C(conductivity, temperature, pressure)
#         assert np.allclose(expected_salinity, result_salinity, rtol=0, atol=1e-4)

#     @pytest.mark.parametrize("cal, cnv, hex, cond_label, use_MV_R", DATASETS)
#     def test_salinity_from_tcp_raw(self, cal, cnv, hex, cond_label, use_MV_R):
#         """Converts data from raw hex values.

#         The unused data loaded from cnv are not needed for the test,
#         but are useful for comparing during a debug session.

#         Args:
#             parameterized from the DATASET const at the top of the file
#         """

#         salinity_expected = cnv.measurements["sal00"].values[800:810]
#         temperature = cnv.measurements["tv290C"].values[800:810]
#         temperature_raw = hex.temperature.values[800:810]
#         temperature_conv = dc.convert_temperature_array(
#             temperature_raw, cal.A0, cal.A1, cal.A2, cal.A3, True, True, use_MV_R
#         )

#         pressure = cnv.measurements["prdM"].values[800:810]
#         pressure_raw = hex.pressure.values[800:810]
#         other_raw = hex.iloc[:, 3].values[800:810]  # either salinity or pressure temperature compensation
#         pressure_conv = dc.convert_pressure_array(
#             pressure_raw,
#             other_raw,
#             True,
#             cal.PA0,
#             cal.PA1,
#             cal.PA2,
#             cal.PTEMPA0,
#             cal.PTEMPA1,
#             cal.PTEMPA2,
#             cal.PTCA0,
#             cal.PTCA1,
#             cal.PTCA2,
#             cal.PTCB0,
#             cal.PTCB1,
#             cal.PTCB2,
#         )

#         conductivity = cnv.measurements[cond_label].values[800:810] * 10
#         conductivity_raw = hex.conductivity.values[800:810]
#         conductivity_conv = (
#             dc.convert_conductivity_array(
#                 conductivity_raw,
#                 temperature_conv,
#                 pressure_conv,
#                 cal.G,
#                 cal.H,
#                 cal.I,
#                 cal.J,
#                 cal.CPCOR,
#                 cal.CTCOR,
#                 cal.WBOTC,
#             )
#             * 10
#         )
#         # can't un-round cnv data, so instead validate that converted hex data rounds to same value
#         salinity_result = np.round(gsw.SP_from_C(conductivity_conv, temperature_conv, pressure_conv), 4)
#         assert np.allclose(salinity_expected, salinity_result, rtol=0, atol=1e-4)


class TestConvertOxygen:
    def test_convert_sbe63_oxygen(self):
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
            ec.oxygen_63_coefs_sn2568,
            ec.thermistor_63_coefs_sn2568,
        )
        assert np.allclose([expected], [result], rtol=0, atol=1e-2)

    def test_convert_sbe43_oxygen(self):
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
            ec.oxygen_43_coefs_sn3287,
            0,
        )
        assert np.allclose([expected], [result], rtol=0, atol=1e-2)

    def test_convert_sbe43_oxygen_from_hex(self):
        # From SBE19plus_01906398_2019_07_15_0033.hex
        raw_oxygen = np.array(
            [
                2.5575,
                2.5586,
                2.5606,
                2.5627,
                2.5638,
                2.5637,
                2.5635,
                2.5629,
                2.5621,
                2.5618,
            ]
        )
        pressure = np.array(
            [-0.012, -0.012, -0.012, -0.012, -0.012, -0.012, -0.012, -0.011, 0.107, 0.351]
        )
        temperature = np.array(
            [
                25.3427,
                25.3408,
                25.3387,
                25.3363,
                25.3341,
                25.3326,
                25.3316,
                25.3302,
                25.3377,
                25.5433,
            ]
        )
        salinity = np.array(
            [0.4373, 0.5592, 0.5865, 0.5095, 0.4621, 0.4119, 0.3936, 0.3463, 4.9297, 6.5098]
        )
        expected = np.array(
            [4.4728, 4.4722, 4.4762, 4.4828, 4.4867, 4.4879, 4.488, 4.488, 4.3707, 4.3148]
        )
        for index in range(len(expected)):
            result = dc._convert_sbe43_oxygen(
                raw_oxygen[index],
                temperature[index],
                pressure[index],
                salinity[index],
                ec.oxygen_43_coefs_sn1686,
                0,
            )
            assert np.allclose([expected[index]], [result], rtol=0, atol=1e-3)

    def test_convert_sbe43_oxygen_from_hex_with_hysteresis(self):
        # TODO: This test is failing. Fix as part of NSI-3061
        # From SBE19plus_01906398_2019_07_15_0033.hex
        # TODO: hysteresis correction only has a real impact on deep data, will need some to better validate this
        raw_oxygen = np.array(
            [
                2.5575,
                2.5586,
                2.5606,
                2.5627,
                2.5638,
                2.5637,
                2.5635,
                2.5629,
                2.5621,
                2.5618,
            ]
        )
        pressure = np.array(
            [-0.012, -0.012, -0.012, -0.012, -0.012, -0.012, -0.012, -0.011, 0.107, 0.351]
        )
        temperature = np.array(
            [
                25.3427,
                25.3408,
                25.3387,
                25.3363,
                25.3341,
                25.3326,
                25.3316,
                25.3302,
                25.3377,
                25.5433,
            ]
        )
        salinity = np.array(
            [0.4373, 0.5592, 0.5865, 0.5095, 0.4621, 0.4119, 0.3936, 0.3463, 4.9297, 6.5098]
        )
        expected = np.array(
            [4.4728, 4.4722, 4.4762, 4.4828, 4.4867, 4.4879, 4.488, 4.488, 4.3707, 4.3148]
        )
        result = dc.convert_sbe43_oxygen(
            raw_oxygen,
            temperature,
            pressure,
            salinity,
            ec.oxygen_43_coefs_sn1686,
            False,
            True,
            1,
            0.25,
        )

        assert np.allclose(result, expected, rtol=0, atol=1e-3)

    def test_convert_sbe43_oxygen_from_hex_with_tau_correction(self):
        # From SBE19plus_01906398_2019_07_15_0033.hex
        raw_oxygen = np.asarray(
            [2.5575, 2.5586, 2.5606, 2.5627, 2.5638, 2.5637, 2.5635, 2.5629, 2.5621, 2.5618]
        )
        pressure = np.asarray(
            [-0.012, -0.012, -0.012, -0.012, -0.012, -0.012, -0.012, -0.011, 0.107, 0.351]
        )
        temperature = np.asarray(
            [
                25.3427,
                25.3408,
                25.3387,
                25.3363,
                25.3341,
                25.3326,
                25.3316,
                25.3302,
                25.3377,
                25.5433,
            ]
        )
        salinity = np.asarray(
            [0.4373, 0.5592, 0.5865, 0.5095, 0.4621, 0.4119, 0.3936, 0.3463, 4.9297, 6.5098]
        )
        expected = [4.4729, 4.4723, 4.4884, 4.4927, 4.4916, 4.4879, 4.4849, 4.4841, 4.3707, 4.3148]

        result = dc.convert_sbe43_oxygen(
            raw_oxygen,
            temperature,
            pressure,
            salinity,
            ec.oxygen_43_coefs_sn1686,
            True,
            False,
            1,
            0.25,
        )

        assert np.allclose(result, expected, rtol=0, atol=1e-4)

    def test_convert_to_mg_per_l(self):
        oxMlPerL = np.array(
            [4.4728, 4.4722, 4.4762, 4.4828, 4.4867, 4.4879, 4.488, 4.488, 4.3707, 4.3148]
        )
        expected = [6.3921, 6.3913, 6.3969, 6.4064, 6.4119, 6.4137, 6.4138, 6.4138, 6.2461, 6.1663]
        result = dc.convert_oxygen_to_mg_per_l(oxMlPerL)
        for index in range(len(expected)):
            assert np.allclose([expected[index]], [result[index]], rtol=0, atol=1e-3)

    def test_convert_to_umol_per_kg(self):
        oxMlPerL = np.array(
            [4.4728, 4.4722, 4.4762, 4.4828, 4.4867, 4.4879, 4.488, 4.488, 4.3707, 4.3148]
        )
        expected = [
            200.3,
            200.254,
            200.427,
            200.735,
            200.916,
            200.979,
            200.984,
            200.991,
            195.064,
            192.356,
        ]
        potentialDensity = np.array(
            [
                -2.7113,
                -2.6188,
                -2.5977,
                -2.6552,
                -2.6903,
                -2.7279,
                -2.7414,
                -2.7768,
                0.6655,
                1.7939,
            ]
        )
        result = dc.convert_oxygen_to_umol_per_kg(oxMlPerL, potentialDensity)
        for index in range(len(expected)):
            assert np.allclose([expected[index]], [result[index]], rtol=0, atol=1e-2)


class TestConvertChlorophylla:
    def test_convert_eco(self):
        rawAnalog = np.array(
            [
                0.0949,
                0.0948,
                0.0960,
                0.0961,
                0.0962,
                0.0959,
                0.1013,
                0.1012,
                0.1015,
                0.1012,
                0.1003,
                0.0999,
                0.0999,
                0.0996,
            ]
        )
        expected = np.array(
            [
                0.2691,
                0.2683,
                0.2798,
                0.2813,
                0.2821,
                0.279,
                0.3332,
                0.3317,
                0.3355,
                0.3317,
                0.3233,
                0.3187,
                0.3195,
                0.3157,
            ]
        )
        for index in range(len(expected)):
            result = dc.convert_eco(rawAnalog[index], ec.chlorophyll_a_coefs_sn6130)
            assert np.allclose([expected[index]], [result], rtol=0, atol=1e-2)


class TestConvertpH:
    def test_convert_sbe18_ph(self):
        rawVolts = np.array(
            [
                2.9507,
                2.9518,
                2.9522,
                2.9517,
                2.9521,
                2.9519,
                2.9517,
                2.9526,
                2.9525,
                2.9537,
                2.953,
                2.9527,
            ]
        )
        temperatureC = np.array(
            [
                26.6927,
                26.6994,
                26.6624,
                26.6122,
                26.5808,
                26.5621,
                26.5503,
                26.5505,
                26.569,
                26.5954,
                26.6182,
                26.6359,
            ]
        )
        expected = np.array(
            [
                8.587,
                8.591,
                8.593,
                8.591,
                8.593,
                8.592,
                8.592,
                8.595,
                8.594,
                8.599,
                8.596,
                8.595,
            ]
        )
        for index in range(len(expected)):
            result = dc.convert_sbe18_ph(rawVolts[index], temperatureC[index], ec.ph_coefs_sn0762)
            assert np.allclose([expected[index]], [result], rtol=0, atol=1e-3)


class TestPARlogarithmic:
    def test_convert_par_logarithmic(self):
        rawVolts = np.array(
            [
                1.176241702,
                1.175249866,
                1.174944686,
                1.176394293,
                0.089951934,
                0.090562294,
                0.090638590,
                0.091172655,
                0.090485999,
                0.090943770,
                0.090409704,
                0.091096360,
                0.093537804,
            ]
        )
        expected = np.array(
            [
                0.81605,
                0.81394,
                0.81330,
                0.81637,
                0.04817,
                0.04824,
                0.04825,
                0.04832,
                0.04823,
                0.04829,
                0.04822,
                0.04831,
                0.04862,
            ]
        )
        for index in range(len(expected)):
            result = dc.convert_par_logarithmic(
                rawVolts[index],
                ec.par_coefs_sn0411,
            )
            assert np.allclose([expected[index]], [result], rtol=0, atol=1e-3)


# class TestContourFromTSP:
#     # Note: this class doesn't actually test anything and is only for debug
#     data_path = './tests/resources/test-data/SBE37SM-RS232_03711722_2015_11_18_subset_derived.asc'
#     data = pd.read_csv(data_path)

#     def test_contour_from_t_s_p_pass(self):
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

    def test_convert_ph_voltage_counts(self):
        internal_ph_volts = dc.convert_ph_voltage_counts(self.internal_ph_counts)
        assert np.allclose(internal_ph_volts, self.expected_internal_ph_volts, atol=1e-6)

    def test_convert_internal_seafet_ph(self):
        internal_ph = dc.convert_internal_seafet_ph(
            ph_counts=self.internal_ph_counts,
            temperature=self.ph_temperature,
            coefs=ec.ph_seafet_internal_coeffs,
        )
        assert np.allclose(internal_ph, self.expected_internal_ph, atol=1e-6)

    def test_convert_external_seafet_ph(self):
        external_ph = dc.convert_external_seafet_ph(
            ph_counts=self.external_ph_counts,
            temperature=self.ph_temperature,
            salinity=35,
            pressure=0,
            coefs=ec.ph_seafet_external_coeffs,
        )
        assert np.allclose(external_ph, self.expected_external_ph, atol=1e-6)
