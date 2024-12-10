"""Data processing unit tests.

"""

# Native imports
from logging import getLogger

# Third-party imports
import numpy as np
import pandas as pd
import pytest

# Sea-Bird imports

# Internal imports
import seabirdscientific.instrument_data as id
import seabirdscientific.processing as dp
import seabirdscientific.conversion as dc
from seabirdscientific.utils import close_enough

logger = getLogger(__name__)

class TestLowPassFilter:
    expected_path = "./tests/resources/test-data/SBE37SM-filtered.asc"
    expected = pd.read_csv(expected_path)["Tv290C"].values

    source_path = "./tests/resources/test-data/SBE37SM-unfiltered.asc"
    source = pd.read_csv(source_path)["Tv290C"].values

    def test_low_pass_filter(self):
        filtered = dp.low_pass_filter(self.source, 10000, sample_interval=120)

        assert np.allclose(self.expected, filtered, rtol=0, atol=1e-4)


class TestAlignCtd:
    expected_data_path = "./tests/resources/test-data/SBE37SM-align.cnv"
    expected_data = id.cnv_to_instrument_data(expected_data_path)
    # expected data was algined with:
    #   sample_interval=120,
    #   tv290C + 12s
    #   cond0S/m + 150s
    #   prdM + 120s
    #   prdE + 0s
    #   sal00 - 240s

    source_data_path = "./tests/resources/test-data/SBE37SM.cnv"
    source_data = id.cnv_to_instrument_data(source_data_path)

    def test_align_ctd_add_simple_pass(self):
        expected_data_path = "./tests/resources/test-data/SBE37SM-align.cnv"
        expected_data = id.cnv_to_instrument_data(expected_data_path)

        source_data_path = "./tests/resources/test-data/SBE37SM.cnv"
        source_data = id.cnv_to_instrument_data(source_data_path).measurements["tv290C"].values

        # Fix last valid sample
        expected_data.measurements["tv290C"].values[len(expected_data.measurements["tv290C"].values) - 1] = -9.99e-29
        result = dp.align_ctd(source_data, 12, 120)
        assert np.allclose(expected_data.measurements["tv290C"].values, result, atol=0.0001)

    def test_align_ctd_add_simple_pass_2(self):
        expected_data_path = "./tests/resources/test-data/SBE37SM-align.cnv"
        expected_data = id.cnv_to_instrument_data(expected_data_path)

        source_data_path = "./tests/resources/test-data/SBE37SM.cnv"
        source_data = id.cnv_to_instrument_data(source_data_path)

        # Fix last valid sample
        expected_data.measurements["cond0S/m"].values[
            len(expected_data.measurements["cond0S/m"].values) - 2
        ] = -9.99e-29
        result = dp.align_ctd(source_data.measurements["cond0S/m"].values, 150, 120)
        assert np.allclose(expected_data.measurements["cond0S/m"].values, result, atol=0.0001)

    def test_align_ctd_add_exact_factor(self):
        expected_data_path = "./tests/resources/test-data/SBE37SM-align.cnv"
        expected_data = id.cnv_to_instrument_data(expected_data_path)

        source_data_path = "./tests/resources/test-data/SBE37SM.cnv"
        source_data = id.cnv_to_instrument_data(source_data_path)

        result = dp.align_ctd(source_data.measurements["prdM"].values, 120, 120)
        assert np.allclose(expected_data.measurements["prdM"].values, result, atol=0.0001)

    def test_align_ctd_no_change(self):
        expected_data_path = "./tests/resources/test-data/SBE37SM-align.cnv"
        expected_data = id.cnv_to_instrument_data(expected_data_path)

        source_data_path = "./tests/resources/test-data/SBE37SM.cnv"
        source_data = id.cnv_to_instrument_data(source_data_path)

        result = dp.align_ctd(source_data.measurements["prdE"].values, 0, 120)
        assert np.allclose(expected_data.measurements["prdE"].values, result, atol=0.0001)

    def test_align_ctd_subtract(self):
        expected_data_path = "./tests/resources/test-data/SBE37SM-align.cnv"
        expected_data = id.cnv_to_instrument_data(expected_data_path)

        source_data_path = "./tests/resources/test-data/SBE37SM.cnv"
        source_data = id.cnv_to_instrument_data(source_data_path)

        # Fix last valid sample
        # expected_data.measurements['sal00'].values[len(expected_data.measurements['sal00'].values) - 2] = -9.99e-29

        result = dp.align_ctd(source_data.measurements["sal00"].values, -240, 120)
        assert np.allclose(expected_data.measurements["sal00"].values, result, atol=0.0001)


class TestCellThermalMass:
    def test_cell_thermal_mass_pass(self):
        expected_data_path = "./tests/resources/test-data/SBE37SM-ctm.cnv"
        expected_data = id.cnv_to_instrument_data(expected_data_path)
        source_data_path = "./tests/resources/test-data/SBE37SM.cnv"
        source_data = id.cnv_to_instrument_data(source_data_path)
        corrected_conductivity = dp.cell_thermal_mass(
            source_data.measurements["tv290C"].values,
            source_data.measurements["cond0S/m"].values,
            0.03,
            7.0,
            120.0,
        )
        assert np.allclose(
            expected_data.measurements["cond0S/m"].values,
            corrected_conductivity,
            atol=0.000001,
        )


class TestLoopEdit:
    def test_loop_edit_pressure_min_velocity_pass(self):
        expected_data = id.cnv_to_instrument_data(
            "./tests/resources/test-data/CAST0002_mod_filt_loop_min_v.cnv"
        )
        data = id.cnv_to_instrument_data("./tests/resources/test-data/CAST0002_mod_filt.cnv")

        dp.loop_edit_pressure(
            pressure=data.measurements["prSM"].values,
            latitude=data.latitude,
            flag=data.measurements["flag"].values,
            sample_interval=data.interval_s,
            min_velocity_type=dp.MinVelocityType.FIXED,
            min_velocity=0.1,
            window_size=3,
            mean_speed_percent=20,
            remove_surface_soak=False,
            min_soak_depth=5,
            max_soak_depth=20,
            use_deck_pressure_offset=False,
            exclude_flags=False,
        )

        expected_flags = expected_data.measurements["flag"].values
        result_flags = data.measurements["flag"].values
        mismatches = sum(expected_flags != result_flags)
        comp = pd.DataFrame({"expected": expected_flags, "result": result_flags})
        error_rows = []

        for n, row in comp.iterrows():
            if row["expected"] != row["result"]:
                error_rows.append(n)

        assert mismatches / len(result_flags) < 0.01

    def test_loop_edit_pressure_min_velocity_remove_soak_pass(self):
        expected_data = id.cnv_to_instrument_data(
            "./tests/resources/test-data/CAST0002_mod_filt_loop_min_v_remove_soak.cnv"
        )
        data = id.cnv_to_instrument_data("./tests/resources/test-data/CAST0002_mod_filt.cnv")

        dp.loop_edit_pressure(
            pressure=data.measurements["prSM"].values,
            latitude=data.latitude,
            flag=data.measurements["flag"].values,
            sample_interval=data.interval_s,
            min_velocity_type=dp.MinVelocityType.FIXED,
            min_velocity=0.1,
            window_size=3,
            mean_speed_percent=20,
            remove_surface_soak=True,
            min_soak_depth=5,
            max_soak_depth=20,
            use_deck_pressure_offset=False,
            exclude_flags=False,
        )

        expected_flags = expected_data.measurements["flag"].values
        result_flags = data.measurements["flag"].values
        mismatches = sum(expected_flags != result_flags)
        comp = pd.DataFrame({"expected": expected_flags, "result": result_flags})
        error_rows = []

        for n, row in comp.iterrows():
            if row["expected"] != row["result"]:
                error_rows.append(n)

        assert mismatches / len(result_flags) < 0.01

    def test_loop_edit_pressure_min_velocity_reset_flags_pass(self):
        expected_data = id.cnv_to_instrument_data(
            "./tests/resources/test-data/CAST0002_mod_filt_loop_min_v.cnv"
        )
        data = id.cnv_to_instrument_data(
            "./tests/resources/test-data/CAST0002_mod_filt_loop_min_v_remove_soak.cnv"
        )

        dp.loop_edit_pressure(
            pressure=data.measurements["prSM"].values,
            latitude=data.latitude,
            flag=data.measurements["flag"].values,
            sample_interval=data.interval_s,
            min_velocity_type=dp.MinVelocityType.FIXED,
            min_velocity=0.1,
            window_size=3,
            mean_speed_percent=20,
            remove_surface_soak=False,
            min_soak_depth=5,
            max_soak_depth=20,
            use_deck_pressure_offset=False,
            exclude_flags=False,
        )

        expected_flags = expected_data.measurements["flag"].values
        result_flags = data.measurements["flag"].values
        mismatches = sum(expected_flags != result_flags)
        comp = pd.DataFrame({"expected": expected_flags, "result": result_flags})
        error_rows = []

        for n, row in comp.iterrows():
            if row["expected"] != row["result"]:
                error_rows.append(n)

        assert mismatches / len(result_flags) < 0.01

    def test_loop_edit_pressure_min_velocity_exclude_flags_pass(self):
        expected_data = id.cnv_to_instrument_data(
            "./tests/resources/test-data/CAST0002_mod_filt_loop_min_v_remove_soak.cnv"
        )
        data = id.cnv_to_instrument_data(
            "./tests/resources/test-data/CAST0002_mod_filt_loop_min_v_exclude_flags_from_remove_soak.cnv"
        )

        dp.loop_edit_pressure(
            pressure=data.measurements["prSM"].values,
            latitude=data.latitude,
            flag=data.measurements["flag"].values,
            sample_interval=data.interval_s,
            min_velocity_type=dp.MinVelocityType.FIXED,
            min_velocity=0.1,
            window_size=3,
            mean_speed_percent=20,
            remove_surface_soak=False,
            min_soak_depth=5,
            max_soak_depth=20,
            use_deck_pressure_offset=False,
            exclude_flags=True,
        )

        expected_flags = expected_data.measurements["flag"].values
        result_flags = data.measurements["flag"].values
        mismatches = sum(expected_flags != result_flags)
        comp = pd.DataFrame({"expected": expected_flags, "result": result_flags})
        error_rows = []

        for n, row in comp.iterrows():
            if row["expected"] != row["result"]:
                error_rows.append(n)

        assert mismatches / len(result_flags) < 0.01

    def test_loop_edit_pressure_mean_speed_percent_remove_soak_pass(self):
        expected_data = id.cnv_to_instrument_data(
            "./tests/resources/test-data/CAST0002_mod_filt_loop_percent_remove_soak.cnv"
        )
        data = id.cnv_to_instrument_data("./tests/resources/test-data/CAST0002_mod_filt.cnv")

        dp.loop_edit_pressure(
            pressure=data.measurements["prSM"].values,
            latitude=data.latitude,
            flag=data.measurements["flag"].values,
            sample_interval=data.interval_s,
            min_velocity_type=dp.MinVelocityType.PERCENT,
            min_velocity=0.1,
            window_size=3,
            mean_speed_percent=20,
            remove_surface_soak=True,
            min_soak_depth=5,
            max_soak_depth=20,
            use_deck_pressure_offset=False,
            exclude_flags=False,
        )

        expected_flags = expected_data.measurements["flag"].values
        result_flags = data.measurements["flag"].values
        mismatches = sum(expected_flags != result_flags)
        comp = pd.DataFrame({"expected": expected_flags, "result": result_flags})
        error_rows = []

        for n, row in comp.iterrows():
            if row["expected"] != row["result"]:
                error_rows.append(n)

        # Less than 2% mismatched flags, outliers appear to conform to spec.
        assert mismatches / len(result_flags) < 0.02

    def test_loop_edit_pressure_mean_speed_percent_pass(self):
        expected_data = id.cnv_to_instrument_data(
            "./tests/resources/test-data/CAST0002_mod_filt_loop_percent.cnv"
        )
        data = id.cnv_to_instrument_data("./tests/resources/test-data/CAST0002_mod_filt.cnv")

        dp.loop_edit_pressure(
            pressure=data.measurements["prSM"].values,
            latitude=data.latitude,
            flag=data.measurements["flag"].values,
            sample_interval=data.interval_s,
            min_velocity_type=dp.MinVelocityType.PERCENT,
            min_velocity=0.1,
            window_size=3,
            mean_speed_percent=20,
            remove_surface_soak=False,
            min_soak_depth=5,
            max_soak_depth=20,
            use_deck_pressure_offset=False,
            exclude_flags=False,
        )

        # SeaSoft is including earliest samples as local maxima during the first
        # downcast which would be discarded when requiring a minimum soak depth
        expected_flags = expected_data.measurements["flag"].values[100:]
        result_flags = data.measurements["flag"].values[100:]
        mismatches = sum(expected_flags != result_flags)
        comp = pd.DataFrame({"expected": expected_flags, "result": result_flags})
        error_rows = []

        for n, row in comp.iterrows():
            if row["expected"] != row["result"]:
                error_rows.append(n)

        # Less than 2% mismatched flags, outliers appear to conform to spec.
        assert mismatches / len(result_flags) < 0.02

    def test_loop_edit_pressure_min_velocity_pass_2(self):
        expected_data = id.cnv_to_instrument_data(
            "./tests/resources/test-data/SBE19plus_loop_edit.cnv"
        )
        data = id.cnv_to_instrument_data(
            "./tests/resources/test-data/SBE19plus.cnv"
        )

        dp.loop_edit_pressure(
            pressure=data.measurements["prdM"].values,
            latitude=data.latitude,
            flag=data.measurements["flag"].values,
            sample_interval=data.interval_s,
            min_velocity_type=dp.MinVelocityType.FIXED,
            min_velocity=0.25,
            window_size=3,
            mean_speed_percent=20,
            remove_surface_soak=True,
            min_soak_depth=5,
            max_soak_depth=20,
            use_deck_pressure_offset=True,
            exclude_flags=True,
        )

        # SeaSoft is including earliest samples as local maxima during the first
        # downcast which would be discarded when requiring a minimum soak depth
        expected_flags = expected_data.measurements["flag"].values[100:]
        result_flags = data.measurements["flag"].values[100:]
        mismatches = sum(expected_flags != result_flags)
        comp = pd.DataFrame({"expected": expected_flags, "result": result_flags})
        error_rows = []

        for n, row in comp.iterrows():
            if row["expected"] != row["result"]:
                error_rows.append(n)

        # Less than 2% mismatched flags, outliers appear to conform to spec.
        assert mismatches / len(result_flags) < 0.02


# class TestBinAverage:


class TestWildEdit:
    def test_wild_edit_pass(self):
        expected_dataset = id.cnv_to_instrument_data(
            "./tests/resources/test-data/19plus_V2_CTD-processing_example_wild_edit.cnv"
        )
        expected_conductivity = expected_dataset.measurements["c0S/m"].values

        dataset = id.cnv_to_instrument_data("./tests/resources/test-data/19plus_V2_CTD-processing_example.cnv")
        conductivity = dataset.measurements["c0S/m"].values
        flags = dataset.measurements["flag"].values
        wild_edit_output = dp.wild_edit(conductivity, flags, 2, 20, 100, 0, False)

        assert np.all(wild_edit_output == expected_conductivity)


# class TestWindowFilter:
#     # file_prefix = "C:\seabirdscientific\documentation\example_data\19plus_V2_CTD-processing_example_window_filter"
#     cnvdata = id.cnv_to_instrument_data(f"{file_prefix}.cnv")
#     pressure = cnvdata.measurements["prdM"].values
#     flags = cnvdata.measurements["flag"].values
#     window_width = 5
#     half_width = 1  # only applies to gaussian
#     offset = 0.25  # only applies to gaussian
#     sample_interval = cnvdata.interval_s  # only applies to gaussian

#     def test_boxcar_filter(self):
#         expected_dataset = id.cnv_to_instrument_data(f"{self.file_prefix}_boxcar_5.cnv")
#         expected_pressure = expected_dataset.measurements["prdM"].values

#         filtered_pressure = dp.window_filter(
#             self.pressure,
#             self.flags,
#             dp.WindowFilterType.BOXCAR,
#             self.window_width,
#             self.cnvdata.interval_s,
#             self.half_width,
#             self.offset,
#             False,
#         )

#         assert close_enough(filtered_pressure, expected_pressure, 3, 1e-12)

#     def test_boxcar_filter_exclude_flags(self):
#         expected_dataset = id.cnv_to_instrument_data(f"{self.file_prefix}_boxcar_5_excluded.cnv")
#         expected_pressure = expected_dataset.measurements["prdM"].values

#         filtered_pressure = dp.window_filter(
#             self.pressure,
#             self.flags,
#             dp.WindowFilterType.BOXCAR,
#             self.window_width,
#             self.cnvdata.interval_s,
#             self.half_width,
#             self.offset,
#             True,
#         )

#         assert close_enough(filtered_pressure, expected_pressure, 3, 1e-12)

#     def test_cosine_filter(self):
#         expected_dataset = id.cnv_to_instrument_data(f"{self.file_prefix}_cosine_5.cnv")
#         expected_pressure = expected_dataset.measurements["prdM"].values

#         filtered_pressure = dp.window_filter(
#             self.pressure,
#             self.flags,
#             dp.WindowFilterType.COSINE,
#             self.window_width,
#             self.cnvdata.interval_s,
#             self.half_width,
#             self.offset,
#             False,
#         )

#         assert close_enough(filtered_pressure, expected_pressure, 3, 1e-12)

#     def test_cosine_filter_exclude_flags(self):
#         expected_dataset = id.cnv_to_instrument_data(f"{self.file_prefix}_cosine_5_excluded.cnv")
#         expected_pressure = expected_dataset.measurements["prdM"].values

#         filtered_pressure = dp.window_filter(
#             self.pressure,
#             self.flags,
#             dp.WindowFilterType.COSINE,
#             self.window_width,
#             self.cnvdata.interval_s,
#             self.half_width,
#             self.offset,
#             True,
#         )

#         assert close_enough(filtered_pressure, expected_pressure, 3, 1e-12)

#     def test_triangle_filter(self):
#         expected_dataset = id.cnv_to_instrument_data(f"{self.file_prefix}_triangle_5.cnv")
#         expected_pressure = expected_dataset.measurements["prdM"].values

#         filtered_pressure = dp.window_filter(
#             self.pressure,
#             self.flags,
#             dp.WindowFilterType.TRIANGLE,
#             self.window_width,
#             self.cnvdata.interval_s,
#             self.half_width,
#             self.offset,
#             False,
#         )

#         assert close_enough(filtered_pressure, expected_pressure, 3, 1e-12)

#     def test_triangle_filter_exclude_flags(self):
#         expected_dataset = id.cnv_to_instrument_data(f"{self.file_prefix}_triangle_5_excluded.cnv")
#         expected_pressure = expected_dataset.measurements["prdM"].values

#         filtered_pressure = dp.window_filter(
#             self.pressure,
#             self.flags,
#             dp.WindowFilterType.TRIANGLE,
#             self.window_width,
#             self.cnvdata.interval_s,
#             self.half_width,
#             self.offset,
#             True,
#         )

#         assert close_enough(filtered_pressure, expected_pressure, 3, 1e-12)

#     def test_gaussian_filter(self):
#         expected_dataset = id.cnv_to_instrument_data(f"{self.file_prefix}_gaussian_5_1_025.cnv")
#         expected_pressure = expected_dataset.measurements["prdM"].values

#         filtered_pressure = dp.window_filter(
#             self.pressure,
#             self.flags,
#             dp.WindowFilterType.GAUSSIAN,
#             self.window_width,
#             self.cnvdata.interval_s,
#             self.half_width,
#             self.offset,
#             False,
#         )

#         assert close_enough(filtered_pressure, expected_pressure, 3, 1e-12)

#     def test_gaussian_filter_exclude_flags(self):
#         expected_dataset = id.cnv_to_instrument_data(f"{self.file_prefix}_gaussian_5_1_0_excluded.cnv")
#         expected_pressure = expected_dataset.measurements["prdM"].values

#         filtered_pressure = dp.window_filter(
#             self.pressure,
#             self.flags,
#             dp.WindowFilterType.GAUSSIAN,
#             self.window_width,
#             self.cnvdata.interval_s,
#             self.half_width,
#             0,
#             True,
#         )

#         assert close_enough(filtered_pressure, expected_pressure, 3, 1e-12)

#     def test_median_filter_exclude_flags(self):
#         # median filter does not match SeaSoft
#         expected_pressure = [
#             7.405,
#             7.722,
#             7.98,
#             8.197,
#             8.286,
#             8.452,
#             8.654,
#             8.779,
#             8.879,
#             8.987,
#             9.169,
#             9.328,
#             9.415,
#             9.605,
#             9.818,
#             10.047,
#             10.305,
#             10.598,
#             10.917,
#             11.082,
#             11.418,
#             11.926,
#             12.265,
#             12.615,
#             12.936,
#             13.401,
#             13.681,
#             13.816,
#             14.051,
#             14.275,
#             14.481,
#             14.663,
#             14.839,
#             15.017,
#             15.192,
#             15.37,
#             15.556,
#             15.751,
#             15.852,
#             16.062,
#             16.408,
#             16.644,
#             16.903,
#             17.18,
#             17.608,
#             17.912,
#             18.068,
#             18.395,
#             18.735,
#             19.08,
#             19.418,
#         ]
#         filtered_pressure = dp.window_filter(
#             self.pressure,
#             self.flags,
#             dp.WindowFilterType.MEDIAN,
#             self.window_width,
#             self.cnvdata.interval_s,
#             self.half_width,
#             0,
#             True,
#         )
#         assert(close_enough(filtered_pressure, expected_pressure, 3, 1e-12))

class TestBuoyancy:
    # Testing data comes from a CalCOFI cruise
    temperature = np.asarray([16.7373, 16.5030, 16.1106, 14.3432, 13.0211, 12.0935, 11.3933, 11.2466, 10.9219, 10.4762, 9.9460])
    pressure = np.asarray([10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110])
    salinity = np.asarray([33.2410, 33.2321, 33.2091, 33.1329, 33.0762, 33.1391, 33.2560, 33.4015, 33.5683, 33.6766, 33.7794])
    expected_N2_win30 = np.asarray([-9.990e-29, 5.7702e-05, 1.9197e-04, 2.6735e-04, 2.1888e-04, 2.1620e-04, 1.7374e-04, 1.5761e-04, 1.6905e-04, 1.6099e-04, -9.990e-29])
    expected_N_win30 = np.asarray([-9.990e-29, 4.35, 7.94, 9.37, 8.48, 8.42, 7.55, 7.19, 7.45, 7.27, -9.990e-29])
    expected_E_win30 = np.asarray([-9.990e-29, 5.8901e-06, 1.9596e-05, 2.7290e-05, 2.2342e-05, 2.2069e-05, 1.7735e-05, 1.6089e-05, 1.7256e-05, 1.6433e-05, -9.990e-29])
    expected_E_pow_8_win30 = np.asarray([-9.990e-29, 589.0, 1959.6, 2729.0, 2234.2, 2206.9, 1773.5, 1608.9, 1725.6, 1643.3, -9.990e-29])

    def test_modern_calc(self):
        # TODO: This test is failing. Fix as part of NSI-3061
        output_dataframe = dp.buoyancy(
            self.temperature,
            self.salinity,
            self.pressure,
            np.asarray([34.034167]), # converted from metadata 34.02.03 N in H,M,S
            np.asarray([121.060556]), # converted from metadata 121 03.38 W in H, M, S
            30, # window size
            True
        )
        print('output for modern calc:')
        print(output_dataframe)

        expected_dataframe = pd.DataFrame()
        expected_dataframe['N2'] = self.expected_N2_win30
        expected_dataframe['N'] = self.expected_N_win30
        expected_dataframe['E'] = self.expected_E_win30
        expected_dataframe['E10^-8'] = self.expected_E_pow_8_win30
        expected_dataframe['N2_diff'] = (output_dataframe.N2 - self.expected_N2_win30)
        expected_dataframe['N2_rpd'] = 100*(output_dataframe.N2 - self.expected_N2_win30)/self.expected_N2_win30

        print('Expected Dataframe')
        print(expected_dataframe)
        print(expected_dataframe.describe())


        # Comparing EOS-80 to TEOS-10 buoyancy calculations.  
        # We do not expect them to agree better than +/-1.5% due to differences in the algorithms
        rel_tol = 0.015  # 1.5%
        assert(output_dataframe['N2'].to_numpy() == pytest.approx(self.expected_N2_win30, rel=rel_tol))
        assert(output_dataframe['N'].to_numpy() == pytest.approx(self.expected_N_win30, rel=rel_tol))
        assert(output_dataframe['E'].to_numpy() == pytest.approx(self.expected_E_win30, rel=rel_tol))
        assert(output_dataframe['E10^-8'].to_numpy() == pytest.approx(self.expected_E_pow_8_win30, rel=rel_tol))

    def test_retro_calc(self): 
        # TODO: This test is failing. Fix as part of NSI-3061
        output_dataframe = dp.buoyancy(
            self.temperature,
            self.salinity,
            self.pressure,
            np.asarray([34.034167]), # converted from metadata 34.02.03 N in H,M,S
            np.asarray([121.060556]), # converted from metadata 121 03.38 W in H, M, S
            30, # window size
            False
        )

        print('output for retro calc:')
        print(output_dataframe)

        expected_dataframe = pd.DataFrame()
        expected_dataframe['N2'] = self.expected_N2_win30
        expected_dataframe['N'] = self.expected_N_win30
        expected_dataframe['E'] = self.expected_E_win30
        expected_dataframe['E10^-8'] = self.expected_E_pow_8_win30
        expected_dataframe['N2_diff'] = (output_dataframe.N2 - self.expected_N2_win30)
        expected_dataframe['N2_rpd'] = 100*(output_dataframe.N2 - self.expected_N2_win30)/self.expected_N2_win30

        print('Expected Dataframe')
        print(expected_dataframe)
        print(expected_dataframe.describe())

        # Comparing SBE Data Processing C++ to local Python results using the same EOS-80 calculations.
        # We expect very very close agreement: << 1% differnce
        rel_tol = 0.0026  # 0.26%
        assert(output_dataframe['N2'].to_numpy() == pytest.approx(self.expected_N2_win30, rel=rel_tol))
        assert(output_dataframe['N'].to_numpy() == pytest.approx(self.expected_N_win30, rel=rel_tol))
        assert(output_dataframe['E'].to_numpy() == pytest.approx(self.expected_E_win30, rel=rel_tol))
        assert(output_dataframe['E10^-8'].to_numpy() == pytest.approx(self.expected_E_pow_8_win30, rel=rel_tol))
