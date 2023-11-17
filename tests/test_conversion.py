"""TODO: test_conversion docstring"""

# Native imports
import importlib.resources

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

test_resources = importlib.resources.files('resources')

DATASETS = [
# cal
# cnv
# hex
# cond_label
# use_MV_R
        (
            cc.SN6130, 
            id.cnv_to_instrument_data(
                test_resources / 'orca-test-data/SBE19plusV2/E8001/E8001.cnv'
            ),
            id.read_hex_file(
                test_resources / 'orca-test-data/SBE19plusV2/E8001/E8001.hex',
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
            True
        ),
        (
            cc.SN03706385, 
            id.cnv_to_instrument_data(
                test_resources / 'orca-test-data/SBE37SM/SBE37SM-6385/SBE37SM-6385.cnv'
            ),
            id.read_hex_file(
                test_resources / 'orca-test-data/SBE37SM/SBE37SM-6385/SBE37SM-6385.hex',
                id.InstrumentType.SBE37SM,
                [
                    id.Sensors.Temperature,
                    id.Sensors.Conductivity,
                    id.Sensors.Pressure,
                ],
            ),
            "cond0S/m",
            False
        ),
        (
            cc.SN03716125, 
            id.cnv_to_instrument_data(
                test_resources / 'orca-test-data/SBE37SM/SBE37SM-RS232_03716125_2017_11_16/SBE37SM-RS232_03716125_2017_11_16.cnv'
            ),
            id.read_hex_file(
                test_resources / 'orca-test-data/SBE37SM/SBE37SM-RS232_03716125_2017_11_16/SBE37SM-RS232_03716125_2017_11_16.hex',
                id.InstrumentType.SBE37SM,
                [
                    id.Sensors.Temperature,
                    id.Sensors.Conductivity,
                    id.Sensors.Pressure,
                ],
            ),
            "cond0S/m",
            False
        ),
    ]


class TestConvertTemperature:
# test temperature raw values
    test_temp_vals = np.array([322798, 322808, 322827, 322838])

    def test_convert_temperature_array_90C(self):
        expected = [20.4459, 20.4451, 20.4436, 20.4427]
        result = dc.convert_temperature_array(
            self.test_temp_vals, cc.SN6130.A0, cc.SN6130.A1, cc.SN6130.A2, cc.SN6130.A3, True, True, True
        )
        assert np.allclose(expected, result, rtol=0, atol=1e-4)

    def test_convert_temperature_array_90CF(self):
        expected = [68.8027, 68.8012, 68.7984, 68.7968]
        result = dc.convert_temperature_array(
            self.test_temp_vals, cc.SN6130.A0, cc.SN6130.A1, cc.SN6130.A2, cc.SN6130.A3, True, False, True
        )
        assert np.allclose(expected, result, rtol=0, atol=1e-4)

    def test_convert_temperature_array_68C(self):
        expected = [20.4508, 20.4500, 20.4485, 20.4476]
        result = dc.convert_temperature_array(
            self.test_temp_vals, cc.SN6130.A0, cc.SN6130.A1, cc.SN6130.A2, cc.SN6130.A3, False, True, True
        )
        assert np.allclose(expected, result, rtol=0, atol=1e-4)
        
    def test_convert_temperature_array_68F(self):
        expected = [68.8115, 68.8100, 68.8073, 68.8057]
        result = dc.convert_temperature_array(
            self.test_temp_vals, cc.SN6130.A0, cc.SN6130.A1, cc.SN6130.A2, cc.SN6130.A3, False, False, True
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
            self.test_pressure_vals, self.test_compensation_vals, False,
            cc.SN6130.PA0, cc.SN6130.PA1, cc.SN6130.PA2, cc.SN6130.PTEMPA0, cc.SN6130.PTEMPA1, cc.SN6130.PTEMPA2, 
            cc.SN6130.PTCA0, cc.SN6130.PTCA1, cc.SN6130.PTCA2, cc.SN6130.PTCB0, cc.SN6130.PTCB1, cc.SN6130.PTCB2
        )
        assert np.allclose(expected, result, rtol=0, atol=1e-3)

    def test_convert_pressure_array_strain_dbar(self):
        expected = [-0.105, -0.106, -0.104, -0.107]
        result = dc.convert_pressure_array(
            self.test_pressure_vals, self.test_compensation_vals, True,
            cc.SN6130.PA0, cc.SN6130.PA1, cc.SN6130.PA2, cc.SN6130.PTEMPA0, cc.SN6130.PTEMPA1, cc.SN6130.PTEMPA2, 
            cc.SN6130.PTCA0, cc.SN6130.PTCA1, cc.SN6130.PTCA2, cc.SN6130.PTCB0, cc.SN6130.PTCB1, cc.SN6130.PTCB2
        )
        assert np.allclose(expected, result, rtol=0, atol=1e-3)


class TestConductivity19plus:
    cnv_path = test_resources / 'orca-test-data/SBE19plusV2/E8001/E8001.cnv'
    expected_data = id.cnv_to_instrument_data(cnv_path)

    hex_path = test_resources / 'orca-test-data/SBE19plusV2/E8001/E8001.hex'
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
        raw['temperature'][0:6], cc.SN6130.A0, cc.SN6130.A1, cc.SN6130.A2, cc.SN6130.A3, False, True, True
    )

    # test_press_data = np.array([-0.105, -0.106, -0.104, -0.107, -0.105, -0.104])
    pressure = dc.convert_pressure_array(
        raw['pressure'][0:6], raw['temperature compensation'][0:6], False,
        cc.SN6130.PA0, cc.SN6130.PA1, cc.SN6130.PA2, cc.SN6130.PTEMPA0, cc.SN6130.PTEMPA1, cc.SN6130.PTEMPA2,
        cc.SN6130.PTCA0, cc.SN6130.PTCA1, cc.SN6130.PTCA2, cc.SN6130.PTCB0, cc.SN6130.PTCB1, cc.SN6130.PTCB2
    )
    
    def test_convert_conductivity_array(self):
        cal = cc.SN6130
        expected = [0.008453, 0.008420, 0.008423, 0.008402, 0.008380, 0.008344]
        result = dc.convert_conductivity_array(
            self.raw['conductivity'][0:6],
            self.temperature,
            self.pressure,
            cal.G, cal.H, cal.I, cal.J, cal.CPCOR, cal.CTCOR, cal.WBOTC
        )
        assert np.allclose(expected, result, rtol=0, atol=1e-6)


    def test_convert_conductivity_val(self):
        cal = cc.SN6130
        expected = [0.008453, 0.008420, 0.008423, 0.008402, 0.008380, 0.008344]
        for index in range(len(expected)):
            result = dc.convert_conductivity_val(
                self.raw['conductivity'][index],
                self.temperature[index],
                self.pressure[index],
                cal.G, cal.H, cal.I, cal.J, cal.CPCOR, cal.CTCOR, cal.WBOTC
            )
            assert np.allclose([expected[index]], [result], rtol=0, atol=1e-6)

class TestConductivity37SM:
    cnv_path = test_resources / 'orca-test-data/SBE37SM/SBE37SM-RS232_03716125_2017_11_16/SBE37SM-RS232_03716125_2017_11_16.cnv'
    expected_data = id.cnv_to_instrument_data(cnv_path)

    hex_path = test_resources / 'orca-test-data/SBE37SM/SBE37SM-RS232_03716125_2017_11_16/SBE37SM-RS232_03716125_2017_11_16.hex'
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
        raw['temperature'][0:6], cc.SN6130.A0, cc.SN6130.A1, cc.SN6130.A2, cc.SN6130.A3, False, True, True
    )

    # test_press_data = np.array([-0.105, -0.106, -0.104, -0.107, -0.105, -0.104])
    pressure = dc.convert_pressure_array(
        raw['pressure'][0:6], raw['temperature compensation'][0:6], False,
        cc.SN03716125.PA0, cc.SN03716125.PA1, cc.SN03716125.PA2, cc.SN03716125.PTEMPA0, cc.SN03716125.PTEMPA1, cc.SN03716125.PTEMPA2,
        cc.SN03716125.PTCA0, cc.SN03716125.PTCA1, cc.SN03716125.PTCA2, cc.SN03716125.PTCB0, cc.SN03716125.PTCB1, cc.SN03716125.PTCB2
    )
    
    def test_convert_conductivity_array(self):
        cal = cc.SN03716125
        expected = [2.711842, 2.715786, 2.715857, 2.715846, 2.715846, 2.715857]
        result = dc.convert_conductivity_array(
            self.raw['conductivity'][0:6],
            self.temperature,
            self.pressure,
            cal.G, cal.H, cal.I, cal.J, cal.CPCOR, cal.CTCOR, cal.WBOTC
        )
        assert np.allclose(expected, result, rtol=0, atol=1e-4)


    def test_convert_conductivity_val(self):
        cal = cc.SN03716125
        expected = [2.711842, 2.715786, 2.715857, 2.715846, 2.715846, 2.715857]
        for index in range(len(expected)):
            result = dc.convert_conductivity_val(
                self.raw['conductivity'][index],
                self.temperature[index],
                self.pressure[index],
                cal.G, cal.H, cal.I, cal.J, cal.CPCOR, cal.CTCOR, cal.WBOTC
            )
            assert np.allclose([expected[index]], [result], rtol=0, atol=1e-4)



class TestDeriveDensity:

    data_path = test_resources / 'orca-test-data/SBE37SM/SBE37SM-RS232_03711722_2015_11_18/SBE37SM-RS232_03711722_2015_11_18_subset_derived.asc'
    data = pd.read_csv(data_path)

    @pytest.mark.parametrize("reference_pressure, expected_column", [
        (0.0, "gsw_sigma0A0"), (1000.0, "gsw_sigma1A0"), (2000.0, "gsw_sigma2A0"),
        (3000.0, "gsw_sigma3A0"), (4000.0, "gsw_sigma4A0")
        ]
    )
    def test_derive_potential_density_from_t_s_p_pass(self, reference_pressure, expected_column):
        temperature_C = self.data['t090C'].values
        salinity_PSU = self.data['sal00'].values
        pressure_dbar = self.data['prM'].values
        expected = self.data[expected_column].values
        result = dc.potential_density_from_t_s_p(
            temperature_C, salinity_PSU, pressure_dbar, reference_pressure=reference_pressure
            )
        # TODO: improve passing condition for all instances of allclose
        assert(np.allclose(result, expected, atol=0.1))


    def test_derive_density_from_t_s_p_pass(self):
        temperature_C = self.data['t090C'].values
        salinity_PSU = self.data['sal00'].values
        pressure_dbar = self.data['prM'].values
        expected = self.data['gsw_densityA0'].values
        result = dc.density_from_t_s_p(temperature_C, salinity_PSU, pressure_dbar)

        assert(np.allclose(result, expected, atol=0.1))


    @pytest.mark.parametrize("reference_pressure, expected_column", [
        (0.0, "gsw_sigma0A0"), (1000.0, "gsw_sigma1A0"), (2000.0, "gsw_sigma2A0"),
        (3000.0, "gsw_sigma3A0"), (4000.0, "gsw_sigma4A0")
        ]
    )
    def test_derive_potential_density_from_t_c_p_pass(self, reference_pressure, expected_column):
        temperature_C = self.data['t090C'].values
        conductivity_mScm = self.data['c0S/m'].values * 10.0
        pressure_dbar = self.data['prM'].values
        expected = self.data[expected_column].values
        result = dc.potential_density_from_t_c_p(
            temperature_C, conductivity_mScm, pressure_dbar, reference_pressure=reference_pressure
            )
        assert(np.allclose(result, expected, atol=0.1))


    def test_derive_density_from_t_c_p_pass(self):
        temperature_C = self.data['t090C'].values
        conductivity_mScm = self.data['c0S/m'].values * 10.0
        pressure_dbar = self.data['prM'].values
        expected = self.data['gsw_densityA0'].values
        result = dc.density_from_t_c_p(temperature_C, conductivity_mScm, pressure_dbar)

        assert(np.allclose(result, expected, atol=0.1))


class TestDepthFromPressure:
    data_path = test_resources / 'orca-test-data/SBE37SM/SBE37SM-6385/SBE37SM-6385.asc'
    data = pd.read_csv(data_path).loc[6:10]
    
    @pytest.mark.parametrize("pressure, pressure_units, depth, depth_units", [
        ("prdM", "dbar", "depSM", "m"),
        ("prdM", "dbar", "depSF", "ft"),
        ("prdE", "psi", "depSM", "m"),
        ("prdE", "psi", "depSF", "ft"),
        ]
    )
    def test_depth_from_pressure_pass(self, pressure, pressure_units, depth, depth_units):
        expected_depth = self.data[depth].values
        pressure = self.data[pressure].values
        result_depth = dc.depth_from_pressure(pressure, 0, depth_units, pressure_units)
        assert(np.allclose(expected_depth, result_depth, atol=0.002))


class TestSalinityFromTCP:

    @pytest.mark.parametrize("cal, cnv, hex, cond_label, use_MV_R", DATASETS)
    def test_salinity_from_tcp(self, cal, cnv, hex, cond_label, use_MV_R):
        expected_salinity = cnv.measurements['sal00'].values[800:810]
        temperature = cnv.measurements['tv290C'].values[800:810]
        pressure = cnv.measurements['prdM'].values[800:810]
        conductivity = cnv.measurements[cond_label].values[800:810] * 10
        result_salinity = gsw.SP_from_C(conductivity, temperature, pressure)
        assert(np.allclose(expected_salinity, result_salinity, rtol=0, atol=1e-4))

    @pytest.mark.parametrize("cal, cnv, hex, cond_label, use_MV_R", DATASETS)
    def test_salinity_from_tcp_raw(self, cal, cnv, hex, cond_label, use_MV_R):
        """Converts data from raw hex values. The unused data loaded from cnv are
        not needed for the test, but are useful for comparing during a debug session.

        Args:
            parameterized from the DATASET const at the top of the file
        """
        salinity_expected = cnv.measurements['sal00'].values[800:810]
        temperature = cnv.measurements['tv290C'].values[800:810]
        temperature_raw = hex.temperature.values[800:810]
        temperature_conv = dc.convert_temperature_array(
            temperature_raw, cal.A0, cal.A1, cal.A2, cal.A3, True, True, use_MV_R
        )

        pressure = cnv.measurements['prdM'].values[800:810]
        pressure_raw = hex.pressure.values[800:810]
        other_raw = hex.iloc[:, 3].values[800:810] # either salinity or pressure temperature compensation
        pressure_conv = dc.convert_pressure_array(
            pressure_raw, other_raw, True,
            cal.PA0, cal.PA1, cal.PA2, cal.PTEMPA0, cal.PTEMPA1, cal.PTEMPA2,
            cal.PTCA0, cal.PTCA1, cal.PTCA2, cal.PTCB0, cal.PTCB1, cal.PTCB2
        )

        conductivity = cnv.measurements[cond_label].values[800:810] * 10
        conductivity_raw = hex.conductivity.values[800:810]
        conductivity_conv = dc.convert_conductivity_array(
            conductivity_raw,
            temperature_conv,
            pressure_conv,
            cal.G, cal.H, cal.I, cal.J, cal.CPCOR, cal.CTCOR, cal.WBOTC
        ) * 10
        # can't un-round cnv data, so instead validate that converted hex data rounds to same value
        salinity_result = np.round(gsw.SP_from_C(conductivity_conv, temperature_conv, pressure_conv), 4)
        assert(np.allclose(salinity_expected, salinity_result, rtol=0, atol=1e-4))

class TestConvertOxygen:
    def test_convert_oxygen(self):
        raw_oxygen = [31.06, 31.66, 32.59, 33.92, 34.82, 35.44]
        pressure = [0, 0, 0, 0, 0, 0]
        temperature = [30, 26, 20, 12, 6, 2]
        salinity = [0, 0, 0, 0, 0, 0]
        expected = [0.706, 0.74, 0.799, 0.892, 1.005, 1.095]
        cal = cc.SN06302568
        for index in range(len(expected)):
            result = dc.convert_oxygen_val(
                raw_oxygen[index],
                temperature[index],
                pressure[index],
                salinity[index],
                cal.a0, cal.a1, cal.a2, cal.b0, cal.b1, cal.c0, cal.c1, cal.c2, cal.e
            )
            assert np.allclose([expected[index]], [result], rtol=0, atol=1e-2)



# class TestContourFromTSP:
#     # Note: this class doesn't actually test anything and is only for debug
#     data_path = test_resources / 'orca-test-data/SBE37SM/SBE37SM-RS232_03711722_2015_11_18/SBE37SM-RS232_03711722_2015_11_18_subset_derived.asc'
#     data = pd.read_csv(data_path)

#     def test_contour_from_t_s_p_pass(self):
#         temperature_C = self.data['t090C'].values
#         salinity_PSU = self.data['sal00'].values
#         pressure_dbar = self.data['prM'].values
#         contour_data = contour.contour_from_t_s_p(
#             temperature_C, salinity_PSU, pressure_dbar
#         )
#         assert(True)
