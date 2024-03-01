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
import sbs.process.cal_coefficients as cc
import sbs.process.conversion as dc
import sbs.process.instrument_data as id

DATASETS = [
    # cal
    # cnv
    # hex
    # cond_label
    # use_MV_R
    (
        cc.SN6130,
        id.cnv_to_instrument_data("./tests/resources/test-data/19plus_V2_CTD-processing_example.cnv"),
        id.read_hex_file("./documentation/example_data/19plus_V2_CTD-processing_example.hex",
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
        ),
        "c0S/m",
        True,
    ),
    # removing until we get a smaller version of this source (NSI-3079)
    # (
    #     cc.SN03706385,
    #     id.cnv_to_instrument_data("./tests/resources/test-data/SBE37SM.cnv"),
    #     id.read_hex_file(
    #         "./tests/resources/test-data/SBE37SM.hex",
    #         id.InstrumentType.SBE37SM,
    #         [
    #             id.Sensors.Temperature,
    #             id.Sensors.Conductivity,
    #             id.Sensors.Pressure,
    #         ],
    #     ),
    #     "cond0S/m",
    #     False,
    # ),
    (
        cc.SN03716125,
        id.cnv_to_instrument_data("./tests/resources/test-data/SBE37SM-RS232_03716125_2017_11_16.cnv"
        ),
        id.read_hex_file("./tests/resources/test-data/SBE37SM-RS232_03716125_2017_11_16.hex",
            id.InstrumentType.SBE37SM,
            [
                id.Sensors.Temperature,
                id.Sensors.Conductivity,
                id.Sensors.Pressure,
            ],
        ),
        "cond0S/m",
        False,
    ),
]


class TestConvertTemperature:
    # test temperature raw values
    test_temp_vals = np.array([322798, 322808, 322827, 322838])

    def test_convert_temperature_array_90C(self):
        expected = [20.4459, 20.4451, 20.4436, 20.4427]
        result = dc.convert_temperature_array(
            self.test_temp_vals,
            cc.SN6130.A0,
            cc.SN6130.A1,
            cc.SN6130.A2,
            cc.SN6130.A3,
            True,
            True,
            True,
        )
        assert np.allclose(expected, result, rtol=0, atol=1e-4)

    def test_convert_temperature_array_90CF(self):
        expected = [68.8027, 68.8012, 68.7984, 68.7968]
        result = dc.convert_temperature_array(
            self.test_temp_vals,
            cc.SN6130.A0,
            cc.SN6130.A1,
            cc.SN6130.A2,
            cc.SN6130.A3,
            True,
            False,
            True,
        )
        assert np.allclose(expected, result, rtol=0, atol=1e-4)

    def test_convert_temperature_array_68C(self):
        expected = [20.4508, 20.4500, 20.4485, 20.4476]
        result = dc.convert_temperature_array(
            self.test_temp_vals,
            cc.SN6130.A0,
            cc.SN6130.A1,
            cc.SN6130.A2,
            cc.SN6130.A3,
            False,
            True,
            True,
        )
        assert np.allclose(expected, result, rtol=0, atol=1e-4)

    def test_convert_temperature_array_68F(self):
        expected = [68.8115, 68.8100, 68.8073, 68.8057]
        result = dc.convert_temperature_array(
            self.test_temp_vals,
            cc.SN6130.A0,
            cc.SN6130.A1,
            cc.SN6130.A2,
            cc.SN6130.A3,
            False,
            False,
            True,
        )
        assert np.allclose(expected, result, rtol=0, atol=1e-4)


class TestConvertPressure:
    # test pressure raw values
    test_pressure_vals = np.array([533539, 533538, 533540, 533537])

    # test pressure temperature compensation raw values
    test_compensation_vals = np.array([20625, 20626, 20626, 20626]) / 13107

    def test_convert_pressure_array_strain_psia(self):
        expected = [-0.153, -0.154, -0.151, -0.156]
        result = dc.convert_pressure_array(
            self.test_pressure_vals,
            self.test_compensation_vals,
            False,
            cc.SN6130.PA0,
            cc.SN6130.PA1,
            cc.SN6130.PA2,
            cc.SN6130.PTEMPA0,
            cc.SN6130.PTEMPA1,
            cc.SN6130.PTEMPA2,
            cc.SN6130.PTCA0,
            cc.SN6130.PTCA1,
            cc.SN6130.PTCA2,
            cc.SN6130.PTCB0,
            cc.SN6130.PTCB1,
            cc.SN6130.PTCB2,
        )
        assert np.allclose(expected, result, rtol=0, atol=1e-3)

    def test_convert_pressure_array_strain_dbar(self):
        expected = [-0.105, -0.106, -0.104, -0.107]
        result = dc.convert_pressure_array(
            self.test_pressure_vals,
            self.test_compensation_vals,
            True,
            cc.SN6130.PA0,
            cc.SN6130.PA1,
            cc.SN6130.PA2,
            cc.SN6130.PTEMPA0,
            cc.SN6130.PTEMPA1,
            cc.SN6130.PTEMPA2,
            cc.SN6130.PTCA0,
            cc.SN6130.PTCA1,
            cc.SN6130.PTCA2,
            cc.SN6130.PTCB0,
            cc.SN6130.PTCB1,
            cc.SN6130.PTCB2,
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
    temperature = dc.convert_temperature_array(
        raw["temperature"][0:6],
        cc.SN6130.A0,
        cc.SN6130.A1,
        cc.SN6130.A2,
        cc.SN6130.A3,
        False,
        True,
        True,
    )

    # test_press_data = np.array([-0.105, -0.106, -0.104, -0.107, -0.105, -0.104])
    pressure = dc.convert_pressure_array(
        raw["pressure"][0:6],
        raw["temperature compensation"][0:6],
        False,
        cc.SN6130.PA0,
        cc.SN6130.PA1,
        cc.SN6130.PA2,
        cc.SN6130.PTEMPA0,
        cc.SN6130.PTEMPA1,
        cc.SN6130.PTEMPA2,
        cc.SN6130.PTCA0,
        cc.SN6130.PTCA1,
        cc.SN6130.PTCA2,
        cc.SN6130.PTCB0,
        cc.SN6130.PTCB1,
        cc.SN6130.PTCB2,
    )

    def test_convert_conductivity_array(self):
        cal = cc.SN6130
        expected = [0.008453, 0.008420, 0.008423, 0.008402, 0.008380, 0.008344]
        result = dc.convert_conductivity_array(
            self.raw["conductivity"][0:6],
            self.temperature,
            self.pressure,
            cal.G,
            cal.H,
            cal.I,
            cal.J,
            cal.CPCOR,
            cal.CTCOR,
            cal.WBOTC,
        )
        assert np.allclose(expected, result, rtol=0, atol=1e-6)

    def test_convert_conductivity_val(self):
        cal = cc.SN6130
        expected = [0.008453, 0.008420, 0.008423, 0.008402, 0.008380, 0.008344]
        for index in range(len(expected)):
            result = dc.convert_conductivity_val(
                self.raw["conductivity"][index],
                self.temperature[index],
                self.pressure[index],
                cal.G,
                cal.H,
                cal.I,
                cal.J,
                cal.CPCOR,
                cal.CTCOR,
                cal.WBOTC,
            )
            assert np.allclose([expected[index]], [result], rtol=0, atol=1e-6)


class TestConductivity37SM:
    cnv_path = ("./tests/resources/test-data/SBE37SM-RS232_03716125_2017_11_16.cnv"
    )
    expected_data = id.cnv_to_instrument_data(cnv_path)

    hex_path = ("./tests/resources/test-data/SBE37SM-RS232_03716125_2017_11_16.hex"
    )
    raw = id.read_hex_file(
        hex_path,
        id.InstrumentType.SBE37SM,
        [
            id.Sensors.Temperature,
            id.Sensors.Conductivity,
            id.Sensors.Pressure,
        ],
    )

    temperature = dc.convert_temperature_array(
        raw["temperature"][0:6],
        cc.SN6130.A0,
        cc.SN6130.A1,
        cc.SN6130.A2,
        cc.SN6130.A3,
        False,
        True,
        True,
    )

    # test_press_data = np.array([-0.105, -0.106, -0.104, -0.107, -0.105, -0.104])
    pressure = dc.convert_pressure_array(
        raw["pressure"][0:6],
        raw["temperature compensation"][0:6],
        False,
        cc.SN03716125.PA0,
        cc.SN03716125.PA1,
        cc.SN03716125.PA2,
        cc.SN03716125.PTEMPA0,
        cc.SN03716125.PTEMPA1,
        cc.SN03716125.PTEMPA2,
        cc.SN03716125.PTCA0,
        cc.SN03716125.PTCA1,
        cc.SN03716125.PTCA2,
        cc.SN03716125.PTCB0,
        cc.SN03716125.PTCB1,
        cc.SN03716125.PTCB2,
    )

    def test_convert_conductivity_array(self):
        cal = cc.SN03716125
        expected = [2.711842, 2.715786, 2.715857, 2.715846, 2.715846, 2.715857]
        result = dc.convert_conductivity_array(
            self.raw["conductivity"][0:6],
            self.temperature,
            self.pressure,
            cal.G,
            cal.H,
            cal.I,
            cal.J,
            cal.CPCOR,
            cal.CTCOR,
            cal.WBOTC,
        )
        assert np.allclose(expected, result, rtol=0, atol=1e-4)

    def test_convert_conductivity_val(self):
        cal = cc.SN03716125
        expected = [2.711842, 2.715786, 2.715857, 2.715846, 2.715846, 2.715857]
        for index in range(len(expected)):
            result = dc.convert_conductivity_val(
                self.raw["conductivity"][index],
                self.temperature[index],
                self.pressure[index],
                cal.G,
                cal.H,
                cal.I,
                cal.J,
                cal.CPCOR,
                cal.CTCOR,
                cal.WBOTC,
            )
            assert np.allclose([expected[index]], [result], rtol=0, atol=1e-4)


class TestDeriveDensity:
    data_path = ("./tests/resources/test-data/SBE37SM-derived.asc"
    )
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
        temperature_C = self.data["t090C"].values
        salinity_PSU = self.data["sal00"].values
        pressure_dbar = self.data["prM"].values
        expected = self.data[expected_column].values
        result = dc.potential_density_from_t_s_p(
            temperature_C,
            salinity_PSU,
            pressure_dbar,
            reference_pressure=reference_pressure,
        )
        # TODO: improve passing condition for all instances of allclose
        assert np.allclose(result, expected, atol=0.1)

    def test_derive_density_from_t_s_p_pass(self):
        temperature_C = self.data["t090C"].values
        salinity_PSU = self.data["sal00"].values
        pressure_dbar = self.data["prM"].values
        expected = self.data["gsw_densityA0"].values
        result = dc.density_from_t_s_p(temperature_C, salinity_PSU, pressure_dbar)

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
        temperature_C = self.data["t090C"].values
        conductivity_mScm = self.data["c0S/m"].values * 10.0
        pressure_dbar = self.data["prM"].values
        expected = self.data[expected_column].values
        result = dc.potential_density_from_t_c_p(
            temperature_C,
            conductivity_mScm,
            pressure_dbar,
            reference_pressure=reference_pressure,
        )
        assert np.allclose(result, expected, atol=0.1)

    def test_derive_density_from_t_c_p_pass(self):
        temperature_C = self.data["t090C"].values
        conductivity_mScm = self.data["c0S/m"].values * 10.0
        pressure_dbar = self.data["prM"].values
        expected = self.data["gsw_densityA0"].values
        result = dc.density_from_t_c_p(temperature_C, conductivity_mScm, pressure_dbar)

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


class TestSalinityFromTCP:
    @pytest.mark.parametrize("cal, cnv, hex, cond_label, use_MV_R", DATASETS)
    def test_salinity_from_tcp(self, cal, cnv, hex, cond_label, use_MV_R):
        expected_salinity = cnv.measurements["sal00"].values[800:810]
        temperature = cnv.measurements["tv290C"].values[800:810]
        pressure = cnv.measurements["prdM"].values[800:810]
        conductivity = cnv.measurements[cond_label].values[800:810] * 10
        result_salinity = gsw.SP_from_C(conductivity, temperature, pressure)
        assert np.allclose(expected_salinity, result_salinity, rtol=0, atol=1e-4)

    @pytest.mark.parametrize("cal, cnv, hex, cond_label, use_MV_R", DATASETS)
    def test_salinity_from_tcp_raw(self, cal, cnv, hex, cond_label, use_MV_R):
        """Converts data from raw hex values. 
        
        The unused data loaded from cnv are not needed for the test, 
        but are useful for comparing during a debug session.

        Args:
            parameterized from the DATASET const at the top of the file
        """
        
        salinity_expected = cnv.measurements["sal00"].values[800:810]
        temperature = cnv.measurements["tv290C"].values[800:810]
        temperature_raw = hex.temperature.values[800:810]
        temperature_conv = dc.convert_temperature_array(
            temperature_raw, cal.A0, cal.A1, cal.A2, cal.A3, True, True, use_MV_R
        )

        pressure = cnv.measurements["prdM"].values[800:810]
        pressure_raw = hex.pressure.values[800:810]
        other_raw = hex.iloc[:, 3].values[800:810]  # either salinity or pressure temperature compensation
        pressure_conv = dc.convert_pressure_array(
            pressure_raw,
            other_raw,
            True,
            cal.PA0,
            cal.PA1,
            cal.PA2,
            cal.PTEMPA0,
            cal.PTEMPA1,
            cal.PTEMPA2,
            cal.PTCA0,
            cal.PTCA1,
            cal.PTCA2,
            cal.PTCB0,
            cal.PTCB1,
            cal.PTCB2,
        )

        conductivity = cnv.measurements[cond_label].values[800:810] * 10
        conductivity_raw = hex.conductivity.values[800:810]
        conductivity_conv = (
            dc.convert_conductivity_array(
                conductivity_raw,
                temperature_conv,
                pressure_conv,
                cal.G,
                cal.H,
                cal.I,
                cal.J,
                cal.CPCOR,
                cal.CTCOR,
                cal.WBOTC,
            )
            * 10
        )
        # can't un-round cnv data, so instead validate that converted hex data rounds to same value
        salinity_result = np.round(gsw.SP_from_C(conductivity_conv, temperature_conv, pressure_conv), 4)
        assert np.allclose(salinity_expected, salinity_result, rtol=0, atol=1e-4)


class TestConvertOxygen:
    def test_convert_sbe63_oxygen(self):
        # without thermistor temperature, this test won't catch the bug fixed in nsi-2656
        raw_oxygen = [31.06, 31.66, 32.59, 33.92, 34.82, 35.44]
        pressure = [0, 0, 0, 0, 0, 0]
        temperature = [30, 26, 20, 12, 6, 2]
        salinity = [0, 0, 0, 0, 0, 0]
        expected = [0.706, 0.74, 0.799, 0.892, 1.005, 1.095]
        cal = cc.SN06302568
        for index in range(len(expected)):
            result = dc.convert_sbe63_oxygen_val(
                raw_oxygen[index],
                temperature[index],
                temperature[index],
                pressure[index],
                salinity[index],
                cal.a0,
                cal.a1,
                cal.a2,
                cal.b0,
                cal.b1,
                cal.c0,
                cal.c1,
                cal.c2,
                cal.e,
                1,
                0,
            )
            assert np.allclose([expected[index]], [result], rtol=0, atol=1e-2)

    def test_convert_sbe43_oxygen(self):
        # From O3287.pdf in the shared calibration folder
        raw_oxygen = [0.725, 0.756, 0.803, 0.874, 0.925, 0.96, 1.332, 1.435, 1.595, 1.81]
        pressure = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        temperature = [2, 6, 12, 20, 26, 30, 2, 6, 12, 20]
        salinity = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        expected = [1.11, 1.12, 1.13, 1.16, 1.18, 1.18, 3.9, 3.9, 3.92, 3.95]
        cal = cc.SN3287

        for index in range(len(expected)):
            result = dc.convert_sbe43_oxygen_val(
                raw_oxygen[index],
                temperature[index],
                pressure[index],
                salinity[index],
                cal.soc,
                cal.v_offset,
                cal.tau_20,
                cal.a,
                cal.b,
                cal.c,
                cal.e,
                cal.d1,
                cal.d2,
                0,
            )
            assert np.allclose([expected[index]], [result], rtol=0, atol=1e-2)

    def test_convert_sbe43_oxygen_from_hex(self):
        # From SBE19plus_01906398_2019_07_15_0033.hex
        raw_oxygen = [2.5575, 2.5586, 2.5606, 2.5627, 2.5638, 2.5637, 2.5635, 2.5629, 2.5621, 2.5618]
        pressure = [
            -0.012,
            -0.012,
            -0.012,
            -0.012,
            -0.012,
            -0.012,
            -0.012,
            -0.011,
            0.107,
            0.351,
        ]
        temperature = [
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
        salinity = [0.4373, 0.5592, 0.5865, 0.5095, 0.4621, 0.4119, 0.3936, 0.3463, 4.9297, 6.5098]
        expected = [4.4728, 4.4722, 4.4762, 4.4828, 4.4867, 4.4879, 4.488, 4.488, 4.3707, 4.3148]
        cal = cc.SN431686
        for index in range(len(expected)):
            result = dc.convert_sbe43_oxygen_val(
                raw_oxygen[index],
                temperature[index],
                pressure[index],
                salinity[index],
                cal.soc,
                cal.v_offset,
                cal.tau_20,
                cal.a,
                cal.b,
                cal.c,
                cal.e,
                cal.d1,
                cal.d2,
                0,
            )
            assert np.allclose([expected[index]], [result], rtol=0, atol=1e-3)

    def test_convert_sbe43_oxygen_from_hex_with_Hysteresis(self):
        # TODO: This test is failing. Fix as part of NSI-3061
        # From SBE19plus_01906398_2019_07_15_0033.hex
        # TODO: hysteresis correction only has a real impact on deep data, will need some to better validate this
        raw_oxygen = [2.5575, 2.5586, 2.5606, 2.5627, 2.5638, 2.5637, 2.5635, 2.5629, 2.5621, 2.5618]
        pressure = [
            -0.012,
            -0.012,
            -0.012,
            -0.012,
            -0.012,
            -0.012,
            -0.012,
            -0.011,
            0.107,
            0.351,
        ]
        temperature = [
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
        salinity = [0.4373, 0.5592, 0.5865, 0.5095, 0.4621, 0.4119, 0.3936, 0.3463, 4.9297, 6.5098]
        expected = [4.4728, 4.4722, 4.4762, 4.4828, 4.4867, 4.4879, 4.488, 4.488, 4.3707, 4.3148]
        cal = cc.SN431686
        result = dc.convert_sbe43_oxygen_array(
            raw_oxygen,
            temperature,
            pressure,
            salinity,
            cal.soc,
            cal.v_offset,
            cal.tau_20,
            cal.a,
            cal.b,
            cal.c,
            cal.e,
            cal.d1,
            cal.d2,
            cal.h1,
            cal.h2,
            cal.h3,
            False,
            True,
            1,
            0.25,
        )

        assert np.allclose(result, expected, rtol=0, atol=1e-3)

    def test_convert_sbe43_oxygen_from_hex_with_tau_correction(self):
        # From SBE19plus_01906398_2019_07_15_0033.hex
        raw_oxygen = np.asarray([2.5575, 2.5586, 2.5606, 2.5627, 2.5638, 2.5637, 2.5635, 2.5629, 2.5621, 2.5618])
        pressure = np.asarray(
            [
                -0.012,
                -0.012,
                -0.012,
                -0.012,
                -0.012,
                -0.012,
                -0.012,
                -0.011,
                0.107,
                0.351,
            ]
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
        salinity = np.asarray([0.4373, 0.5592, 0.5865, 0.5095, 0.4621, 0.4119, 0.3936, 0.3463, 4.9297, 6.5098])
        expected = [4.4729, 4.4723, 4.4884, 4.4927, 4.4916, 4.4879, 4.4849, 4.4841, 4.3707, 4.3148]
        cal = cc.SN431686

        result = dc.convert_sbe43_oxygen_array(
            raw_oxygen,
            temperature,
            pressure,
            salinity,
            cal.soc,
            cal.v_offset,
            cal.tau_20,
            cal.a,
            cal.b,
            cal.c,
            cal.e,
            cal.d1,
            cal.d2,
            cal.h1,
            cal.h2,
            cal.h3,
            True,
            False,
            1,
            0.25,
        )

        assert np.allclose(result, expected, rtol=0, atol=1e-4)

    def test_convert_to_mg_per_l(self):
        oxMlPerL = np.array(
            [
                4.4728,
                4.4722,
                4.4762,
                4.4828,
                4.4867,
                4.4879,
                4.488,
                4.488,
                4.3707,
                4.3148,
            ]
        )
        expected = [
            6.3921,
            6.3913,
            6.3969,
            6.4064,
            6.4119,
            6.4137,
            6.4138,
            6.4138,
            6.2461,
            6.1663,
        ]
        result = dc.convert_oxygen_to_mg_per_l(oxMlPerL)
        for index in range(len(expected)):
            assert np.allclose([expected[index]], [result[index]], rtol=0, atol=1e-3)

    def test_convert_to_umol_per_kg(self):
        oxMlPerL = np.array(
            [
                4.4728,
                4.4722,
                4.4762,
                4.4828,
                4.4867,
                4.4879,
                4.488,
                4.488,
                4.3707,
                4.3148,
            ]
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
    def test_convert_ECO_chlorophylla(self):
        rawAnalog = [
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
        expected = [
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
        cal = cc.SN6130
        for index in range(len(expected)):
            result = dc.convert_ECO_chlorophylla_val(rawAnalog[index], cal.ScaleFactorChla, cal.Vblank)
            assert np.allclose([expected[index]], [result], rtol=0, atol=1e-2)


class TestConvertTurbidity:
    def test_convert_ECO_turbidity(self):
        rawAnalog = [
            0.0787,
            0.079,
            0.0831,
            0.0829,
            0.0835,
            0.0833,
            0.0825,
            0.082,
            0.082,
            0.0822,
            0.0812,
            0.0806,
            0.0813,
            0.0816,
        ]
        expected = [
            0.0983,
            0.1002,
            0.1204,
            0.1197,
            0.1223,
            0.1216,
            0.1174,
            0.1151,
            0.1151,
            0.1158,
            0.1109,
            0.1082,
            0.1117,
            0.1132,
        ]
        cal = cc.SN6130
        for index in range(len(expected)):
            result = dc.convert_ECO_turbidity_val(rawAnalog[index], cal.ScaleFactorTurbidity, cal.DarkVoltage)
            assert np.allclose([expected[index]], [result], rtol=0, atol=1e-2)


# class TestContourFromTSP:
#     # Note: this class doesn't actually test anything and is only for debug
#     data_path = './tests/resources/test-data/SBE37SM-RS232_03711722_2015_11_18_subset_derived.asc'
#     data = pd.read_csv(data_path)

#     def test_contour_from_t_s_p_pass(self):
#         temperature_C = self.data['t090C'].values
#         salinity_PSU = self.data['sal00'].values
#         pressure_dbar = self.data['prM'].values
#         contour_data = contour.contour_from_t_s_p(
#             temperature_C, salinity_PSU, pressure_dbar
#         )
#         assert(True)
