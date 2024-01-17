"""TODO: test_processing docstring"""

# Native imports
import importlib.resources
from logging import getLogger

# Third-party imports
import numpy as np
import pandas as pd

# Sea-Bird imports

# Internal imports
from sbs.process import instrument_data as id
from sbs.process import processing as dp
from sbs.process import conversion as dc
from sbs.process.utils import close_enough

logger = getLogger(__name__)

test_resources = importlib.resources.files("resources")
microcat6385_resources = test_resources / "orca-test-data" / "SBE37SM" / "SBE37SM-6385"


class TestLowPassFilter:
    expected_path = microcat6385_resources / "SBE37SM-6385_filtered_10000.asc"
    expected = pd.read_csv(expected_path)["Tv290C"].values

    source_path = microcat6385_resources / "SBE37SM-6385_unfiltered.asc"
    source = pd.read_csv(source_path)["Tv290C"].values

    def test_low_pass_filter(self):
        filtered = dp.low_pass_filter(self.source, 10000, sample_interval=120)

        assert np.allclose(self.expected, filtered, rtol=0, atol=1e-4)


class TestAlignCtd:
    expected_data_path = microcat6385_resources / "SBE37SM-6385_align.cnv"
    expected_data = id.cnv_to_instrument_data(expected_data_path)
    # expected data was algined with:
    #   sample_interval=120,
    #   tv290C + 12s
    #   cond0S/m + 150s
    #   prdM + 120s
    #   prdE + 0s
    #   sal00 - 240s

    source_data_path = microcat6385_resources / "SBE37SM-6385.cnv"
    source_data = id.cnv_to_instrument_data(source_data_path)

    def test_align_ctd_add_simple_pass(self):
        expected_data_path = microcat6385_resources / "SBE37SM-6385_align.cnv"
        expected_data = id.cnv_to_instrument_data(expected_data_path)

        source_data_path = microcat6385_resources / "SBE37SM-6385.cnv"
        source_data = id.cnv_to_instrument_data(source_data_path).measurements["tv290C"].values

        # Fix last valid sample
        expected_data.measurements["tv290C"].values[len(expected_data.measurements["tv290C"].values) - 1] = -9.99e-29
        result = dp.align_ctd(source_data, 12, 120)
        assert np.allclose(expected_data.measurements["tv290C"].values, result, atol=0.0001)

    def test_align_ctd_add_simple_pass_2(self):
        expected_data_path = microcat6385_resources / "SBE37SM-6385_align.cnv"
        expected_data = id.cnv_to_instrument_data(expected_data_path)

        source_data_path = microcat6385_resources / "SBE37SM-6385.cnv"
        source_data = id.cnv_to_instrument_data(source_data_path)

        # Fix last valid sample
        expected_data.measurements["cond0S/m"].values[
            len(expected_data.measurements["cond0S/m"].values) - 2
        ] = -9.99e-29
        result = dp.align_ctd(source_data.measurements["cond0S/m"].values, 150, 120)
        assert np.allclose(expected_data.measurements["cond0S/m"].values, result, atol=0.0001)

    def test_align_ctd_add_exact_factor(self):
        expected_data_path = microcat6385_resources / "SBE37SM-6385_align.cnv"
        expected_data = id.cnv_to_instrument_data(expected_data_path)

        source_data_path = microcat6385_resources / "SBE37SM-6385.cnv"
        source_data = id.cnv_to_instrument_data(source_data_path)

        result = dp.align_ctd(source_data.measurements["prdM"].values, 120, 120)
        assert np.allclose(expected_data.measurements["prdM"].values, result, atol=0.0001)

    def test_align_ctd_no_change(self):
        expected_data_path = microcat6385_resources / "SBE37SM-6385_align.cnv"
        expected_data = id.cnv_to_instrument_data(expected_data_path)

        source_data_path = microcat6385_resources / "SBE37SM-6385.cnv"
        source_data = id.cnv_to_instrument_data(source_data_path)

        result = dp.align_ctd(source_data.measurements["prdE"].values, 0, 120)
        assert np.allclose(expected_data.measurements["prdE"].values, result, atol=0.0001)

    def test_align_ctd_subtract(self):
        expected_data_path = microcat6385_resources / "SBE37SM-6385_align.cnv"
        expected_data = id.cnv_to_instrument_data(expected_data_path)

        source_data_path = microcat6385_resources / "SBE37SM-6385.cnv"
        source_data = id.cnv_to_instrument_data(source_data_path)

        # Fix last valid sample
        # expected_data.measurements['sal00'].values[len(expected_data.measurements['sal00'].values) - 2] = -9.99e-29

        result = dp.align_ctd(source_data.measurements["sal00"].values, -240, 120)
        assert np.allclose(expected_data.measurements["sal00"].values, result, atol=0.0001)


class TestCellThermalMass:
    def test_cell_thermal_mass_pass(self):
        expected_data_path = microcat6385_resources / "SBE37SM-6385_ctm.cnv"
        expected_data = id.cnv_to_instrument_data(expected_data_path)
        source_data_path = microcat6385_resources / "SBE37SM-6385.cnv"
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
            test_resources / "orca-test-data/SBE19/CAST0002/CAST0002_mod_filt_loop_min_v.cnv"
        )
        data = id.cnv_to_instrument_data(test_resources / "orca-test-data/SBE19/CAST0002/CAST0002_mod_filt.cnv")

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
            test_resources / "orca-test-data/SBE19/CAST0002/CAST0002_mod_filt_loop_min_v_remove_soak.cnv"
        )
        data = id.cnv_to_instrument_data(test_resources / "orca-test-data/SBE19/CAST0002/CAST0002_mod_filt.cnv")

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
            test_resources / "orca-test-data/SBE19/CAST0002/CAST0002_mod_filt_loop_min_v.cnv"
        )
        data = id.cnv_to_instrument_data(
            test_resources / "orca-test-data/SBE19/CAST0002/CAST0002_mod_filt_loop_min_v_remove_soak.cnv"
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
            test_resources / "orca-test-data/SBE19/CAST0002/CAST0002_mod_filt_loop_min_v_remove_soak.cnv"
        )
        data = id.cnv_to_instrument_data(
            test_resources
            / "orca-test-data/SBE19/CAST0002/CAST0002_mod_filt_loop_min_v_exclude_flags_from_remove_soak.cnv"
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
            test_resources / "orca-test-data/SBE19/CAST0002/CAST0002_mod_filt_loop_percent_remove_soak.cnv"
        )
        data = id.cnv_to_instrument_data(test_resources / "orca-test-data/SBE19/CAST0002/CAST0002_mod_filt.cnv")

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
            test_resources / "orca-test-data/SBE19/CAST0002/CAST0002_mod_filt_loop_percent.cnv"
        )
        data = id.cnv_to_instrument_data(test_resources / "orca-test-data/SBE19/CAST0002/CAST0002_mod_filt.cnv")

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
            test_resources
            / "orca-test-data/SBE19plusV2/SBE19plus_01906398_2019_07_15/SBE19plus_01906398_2019_07_15_0033_loop_edit.cnv"
        )
        data = id.cnv_to_instrument_data(
            test_resources
            / "orca-test-data/SBE19plusV2/SBE19plus_01906398_2019_07_15/SBE19plus_01906398_2019_07_15_0033.cnv"
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
            test_resources / r"orca-test-data\SBE19plusV2\E8001\E8001_wild_edit.cnv"
        )
        expected_output = expected_dataset.measurements["prdM"].values

        dataset = id.cnv_to_instrument_data(test_resources / r"orca-test-data\SBE19plusV2\E8001\E8001.cnv")
        pressure = dataset.measurements["prdM"].values
        flags = dataset.measurements["flag"].values

        wild_edit_output = dp.wild_edit(pressure, flags, 1, 3, 100, 0, False)

        assert np.all(wild_edit_output == expected_output)


class TestWindowFilter:
    file_prefix = test_resources / "orca-test-data/SBE19plusV2/E8001/E8001_window_filter"
    cnvdata = id.cnv_to_instrument_data(f"{file_prefix}.cnv")
    pressure = cnvdata.measurements["prdM"].values
    flags = cnvdata.measurements["flag"].values
    window_width = 5
    half_width = 1  # only applies to gaussian
    offset = 0.25  # only applies to gaussian
    sample_interval = cnvdata.interval_s  # only applies to gaussian

    def test_boxcar_filter(self):
        expected_dataset = id.cnv_to_instrument_data(f"{self.file_prefix}_boxcar_5.cnv")
        expected_pressure = expected_dataset.measurements["prdM"].values

        filtered_pressure = dp.window_filter(
            self.pressure,
            self.flags,
            dp.WindowFilterType.BOXCAR,
            self.window_width,
            self.cnvdata.interval_s,
            self.half_width,
            self.offset,
            False,
        )

        assert close_enough(filtered_pressure, expected_pressure, 3, 1e-12)

    def test_boxcar_filter_exclude_flags(self):
        expected_dataset = id.cnv_to_instrument_data(f"{self.file_prefix}_boxcar_5_excluded.cnv")
        expected_pressure = expected_dataset.measurements["prdM"].values

        filtered_pressure = dp.window_filter(
            self.pressure,
            self.flags,
            dp.WindowFilterType.BOXCAR,
            self.window_width,
            self.cnvdata.interval_s,
            self.half_width,
            self.offset,
            True,
        )

        assert close_enough(filtered_pressure, expected_pressure, 3, 1e-12)

    def test_cosine_filter(self):
        expected_dataset = id.cnv_to_instrument_data(f"{self.file_prefix}_cosine_5.cnv")
        expected_pressure = expected_dataset.measurements["prdM"].values

        filtered_pressure = dp.window_filter(
            self.pressure,
            self.flags,
            dp.WindowFilterType.COSINE,
            self.window_width,
            self.cnvdata.interval_s,
            self.half_width,
            self.offset,
            False,
        )

        assert close_enough(filtered_pressure, expected_pressure, 3, 1e-12)

    def test_cosine_filter_exclude_flags(self):
        expected_dataset = id.cnv_to_instrument_data(f"{self.file_prefix}_cosine_5_excluded.cnv")
        expected_pressure = expected_dataset.measurements["prdM"].values

        filtered_pressure = dp.window_filter(
            self.pressure,
            self.flags,
            dp.WindowFilterType.COSINE,
            self.window_width,
            self.cnvdata.interval_s,
            self.half_width,
            self.offset,
            True,
        )

        assert close_enough(filtered_pressure, expected_pressure, 3, 1e-12)

    def test_triangle_filter(self):
        expected_dataset = id.cnv_to_instrument_data(f"{self.file_prefix}_triangle_5.cnv")
        expected_pressure = expected_dataset.measurements["prdM"].values

        filtered_pressure = dp.window_filter(
            self.pressure,
            self.flags,
            dp.WindowFilterType.TRIANGLE,
            self.window_width,
            self.cnvdata.interval_s,
            self.half_width,
            self.offset,
            False,
        )

        assert close_enough(filtered_pressure, expected_pressure, 3, 1e-12)

    def test_triangle_filter_exclude_flags(self):
        expected_dataset = id.cnv_to_instrument_data(f"{self.file_prefix}_triangle_5_excluded.cnv")
        expected_pressure = expected_dataset.measurements["prdM"].values

        filtered_pressure = dp.window_filter(
            self.pressure,
            self.flags,
            dp.WindowFilterType.TRIANGLE,
            self.window_width,
            self.cnvdata.interval_s,
            self.half_width,
            self.offset,
            True,
        )

        assert close_enough(filtered_pressure, expected_pressure, 3, 1e-12)

    def test_gaussian_filter(self):
        expected_dataset = id.cnv_to_instrument_data(f"{self.file_prefix}_gaussian_5_1_025.cnv")
        expected_pressure = expected_dataset.measurements["prdM"].values

        filtered_pressure = dp.window_filter(
            self.pressure,
            self.flags,
            dp.WindowFilterType.GAUSSIAN,
            self.window_width,
            self.cnvdata.interval_s,
            self.half_width,
            self.offset,
            False,
        )

        assert close_enough(filtered_pressure, expected_pressure, 3, 1e-12)

    def test_gaussian_filter_exclude_flags(self):
        expected_dataset = id.cnv_to_instrument_data(f"{self.file_prefix}_gaussian_5_1_0_excluded.cnv")
        expected_pressure = expected_dataset.measurements["prdM"].values

        filtered_pressure = dp.window_filter(
            self.pressure,
            self.flags,
            dp.WindowFilterType.GAUSSIAN,
            self.window_width,
            self.cnvdata.interval_s,
            self.half_width,
            0,
            True,
        )

        assert close_enough(filtered_pressure, expected_pressure, 3, 1e-12)

    def test_median_filter_exclude_flags(self):
        # median filter does not match SeaSoft
        expected_pressure = [
            7.405,
            7.722,
            7.98,
            8.197,
            8.286,
            8.452,
            8.654,
            8.779,
            8.879,
            8.987,
            9.169,
            9.328,
            9.415,
            9.605,
            9.818,
            10.047,
            10.305,
            10.598,
            10.917,
            11.082,
            11.418,
            11.926,
            12.265,
            12.615,
            12.936,
            13.401,
            13.681,
            13.816,
            14.051,
            14.275,
            14.481,
            14.663,
            14.839,
            15.017,
            15.192,
            15.37,
            15.556,
            15.751,
            15.852,
            16.062,
            16.408,
            16.644,
            16.903,
            17.18,
            17.608,
            17.912,
            18.068,
            18.395,
            18.735,
            19.08,
            19.418,
        ]
        filtered_pressure = dp.window_filter(
            self.pressure,
            self.flags,
            dp.WindowFilterType.MEDIAN,
            self.window_width,
            self.cnvdata.interval_s,
            self.half_width,
            0,
            True,
        )
        assert close_enough(filtered_pressure, expected_pressure, 3, 1e-12)
