"""Data processing unit tests."""

# Native imports
from logging import getLogger
from pathlib import Path

# Third-party imports
import numpy as np
import pandas as pd
import pytest

# Sea-Bird imports

# Internal imports
import seabirdscientific.conversion as sc
import seabirdscientific.instrument_data as si
import seabirdscientific.processing as sp
from seabirdscientific.utils import close_enough, get_tolerance

test_data = Path("./tests/resources/test-data")
logger = getLogger(__name__)


class TestLowPassFilter:
    def test_low_pass_filter(self, request):
        expected_path = test_data / "SBE37SM-filtered.asc"
        expected = pd.read_csv(expected_path)["Tv290C"].values

        source_path = test_data / "SBE37SM-unfiltered.asc"
        source = pd.read_csv(source_path)["Tv290C"].values

        filtered = sp.low_pass_filter(source, 10000, sample_interval=120)

        request.node.return_value = filtered.tolist()

        assert np.allclose(expected, filtered, rtol=0, atol=1e-4)


class TestAlignCtdFathomCases:
    # Test align CTD with same test data as Fathom
    def test_align_ctd_basic_forward(self, request):
        badFlagValue = -9.99e-29

        align_cases = [
            # Basic forward align cases
            [
                np.array([0, 15, 30, 45, 60, 75, 90, 105]),  # original
                15,  # offset
                15,  # sample interval
                np.array([15, 30, 45, 60, 75, 90, 105, badFlagValue]),  # expected
            ],
            [
                np.array([0, 15, 30, 45, 60, 75, 90, 105]),
                10,
                15,
                np.array([10, 25, 40, 55, 70, 85, 100, badFlagValue]),
            ],
            [
                np.array([0, 15, 30, 45, 60, 75, 90, 105]),
                5,
                15,
                np.array([5, 20, 35, 50, 65, 80, 95, badFlagValue]),
            ],
            [
                np.array([0, 15, 30, 45, 60, 75, 90, 105]),
                0,
                15,
                np.array([0, 15, 30, 45, 60, 75, 90, 105]),
            ],
            [
                np.array([0, 15, 30, 45, 60, 75, 90, 105]),
                11.575,
                15,
                np.array([11.575, 26.575, 41.575, 56.575, 71.575, 86.575, 101.575, badFlagValue]),
            ],
            # Advanced forward align cases
            [
                np.array([0, 15, 30, 45, 60, 75, 90, 105]),
                30,
                15,
                np.array([30, 45, 60, 75, 90, 105, badFlagValue, badFlagValue]),
            ],
            [
                np.array([0, 15, 30, 45, 60, 75, 90, 105]),
                20,
                15,
                np.array([20, 35, 50, 65, 80, 95, badFlagValue, badFlagValue]),
            ],
            [
                np.array([0, 15, 30, 45, 60, 75, 90, 105]),
                106,
                15,
                np.array([badFlagValue] * 8),
            ],
            [
                np.array([0, 15, 30, 45, 60, 75, 90, 105]),
                60,
                15,
                np.array(
                    [60, 75, 90, 105, badFlagValue, badFlagValue, badFlagValue, badFlagValue]
                ),
            ],
            [
                np.array([0, 15, 30, 45, 60, 75, 90, 105]),
                130,
                15,
                np.array([badFlagValue] * 8),
            ],
            # Basic backwards align cases
            [
                np.array([0, 15, 30, 45, 60, 75, 90, 105]),
                -15,
                15,
                np.array([badFlagValue, 0, 15, 30, 45, 60, 75, 90]),
            ],
            [
                np.array([0, 15, 30, 45, 60, 75, 90, 105]),
                -10,
                15,
                np.array([badFlagValue, 5, 20, 35, 50, 65, 80, 95]),
            ],
            [
                np.array([0, 15, 30, 45, 60, 75, 90, 105]),
                -5,
                15,
                np.array([badFlagValue, 10, 25, 40, 55, 70, 85, 100]),
            ],
            # Advanced backwards align cases
            [
                np.array([0, 15, 30, 45, 60, 75, 90, 105]),
                -30,
                15,
                np.array([badFlagValue, badFlagValue, 0, 15, 30, 45, 60, 75]),
            ],
            [
                np.array([0, 15, 30, 45, 60, 75, 90, 105]),
                -20,
                15,
                np.array([badFlagValue, badFlagValue, 10, 25, 40, 55, 70, 85]),
            ],
            [
                np.array([0, 15, 30, 45, 60, 75, 90, 105]),
                -106,
                15,
                np.array([badFlagValue] * 8),
            ],
            [
                np.array([0, 15, 30, 45, 60, 75, 90, 105]),
                -130,
                15,
                np.array([badFlagValue] * 8),
            ],
        ]

        for case in align_cases:
            result = sp.align_ctd(case[0], case[1], case[2])
            assert np.array_equal(case[3], result)


class TestAlignCtd:
    expected_data_path = test_data / "SBE37SM-align.cnv"
    # expected_data = si.read_cnv_file(expected_data_path)
    # expected data was algined with:
    #   sample_interval=120,
    #   tv290C + 12s
    #   cond0S/m + 150s
    #   prdM + 120s
    #   prdE + 0s
    #   sal00 - 240s

    source_data_path = test_data / "SBE37SM.cnv"
    # source_data = si.read_cnv_file(source_data_path)

    @pytest.fixture
    def expected_data(self):
        return si.read_cnv_file(self.expected_data_path)

    @pytest.fixture
    def source_data(self):
        return si.read_cnv_file(self.source_data_path)

    def test_align_ctd_add_simple_pass(self, expected_data, source_data, request):
        # Fix last valid sample
        expected_data["tv290C"].values[len(expected_data["tv290C"].values) - 1] = -9.99e-29
        result = sp.align_ctd(source_data["tv290C"].values, 12, 120)
        request.node.return_value = result.tolist()
        assert np.allclose(expected_data["tv290C"].values, result, atol=0.0001)

    def test_align_ctd_add_simple_pass_2(self, expected_data, source_data, request):
        # Fix last valid sample
        expected_data["cond0S/m"].values[len(expected_data["cond0S/m"].values) - 2] = -9.99e-29
        result = sp.align_ctd(source_data["cond0S/m"].values, 150, 120)
        request.node.return_value = result.tolist()
        assert np.allclose(expected_data["cond0S/m"].values, result, atol=0.0001)

    def test_align_ctd_add_exact_factor(self, expected_data, source_data, request):
        result = sp.align_ctd(source_data["prdM"].values, 120, 120)
        request.node.return_value = result.tolist()
        assert np.allclose(expected_data["prdM"].values, result, atol=0.0001)

    def test_align_ctd_no_change(self, expected_data, source_data, request):
        result = sp.align_ctd(source_data["prdE"].values, 0, 120)
        request.node.return_value = result.tolist()
        assert np.allclose(expected_data["prdE"].values, result, atol=0.0001)

    def test_align_ctd_subtract(self, expected_data, source_data, request):
        result = sp.align_ctd(source_data["sal00"].values, -240, 120)
        request.node.return_value = result.tolist()
        assert np.allclose(expected_data["sal00"].values, result, atol=0.0001)


class TestCellThermalMass:
    def test_cell_thermal_mass_pass(self, request):
        expected_data_path = test_data / "SBE37SM-ctm.cnv"
        expected_data = si.read_cnv_file(expected_data_path)
        source_data_path = test_data / "SBE37SM.cnv"
        source_data = si.read_cnv_file(source_data_path)
        corrected_conductivity = sp.cell_thermal_mass(
            source_data["tv290C"].values,
            source_data["cond0S/m"].values,
            0.03,
            7.0,
            120.0,
        )
        request.node.return_value = corrected_conductivity.tolist()
        assert np.allclose(
            expected_data["cond0S/m"].values,
            corrected_conductivity,
            atol=0.000001,
        )


class TestLoopEdit:
    def test_loop_edit_pressure_min_velocity_pass(self):
        expected_data = si.read_cnv_file(test_data / "CAST0002_mod_filt_loop_min_v.cnv")
        data = si.read_cnv_file(test_data / "CAST0002_mod_filt.cnv")

        result_flags = sp.loop_edit_pressure(
            pressure=data["prSM"].values,
            latitude=0,
            flag=data["flag"].values,
            sample_interval=data.attrs["sample_interval"],
            min_velocity_type="fixed",
            min_velocity=0.1,
            window_size=3,
            mean_speed_percent=20,
            remove_surface_soak=False,
            min_soak_depth=5,
            max_soak_depth=20,
            use_deck_pressure_offset=False,
            exclude_flags=False,
        )

        expected_flags = expected_data["flag"].values
        assert np.all(expected_flags == result_flags)

    def test_loop_edit_pressure_min_velocity_remove_soak_pass(self, request):
        expected_data = si.read_cnv_file(
            test_data / "CAST0002_mod_filt_loop_min_v_remove_soak.cnv"
        )
        data = si.read_cnv_file(test_data / "CAST0002_mod_filt.cnv")

        result_flags = sp.loop_edit_pressure(
            pressure=data["prSM"].values,
            latitude=0,
            flag=data["flag"].values,
            sample_interval=data.attrs["sample_interval"],
            min_velocity_type="fixed",
            min_velocity=0.1,
            window_size=3,
            mean_speed_percent=20,
            remove_surface_soak=True,
            min_soak_depth=5,
            max_soak_depth=20,
            use_deck_pressure_offset=False,
            exclude_flags=False,
        )

        expected_flags = expected_data["flag"].values
        assert np.all(expected_flags == result_flags)

    def test_loop_edit_pressure_min_velocity_reset_flags_pass(self, request):
        expected_data = si.read_cnv_file(test_data / "CAST0002_mod_filt_loop_min_v.cnv")
        data = si.read_cnv_file(test_data / "CAST0002_mod_filt_loop_min_v_remove_soak.cnv")

        result_flags = sp.loop_edit_pressure(
            pressure=data["prSM"].values,
            latitude=0,
            flag=data["flag"].values,
            sample_interval=data.attrs["sample_interval"],
            min_velocity_type="fixed",
            min_velocity=0.1,
            window_size=3,
            mean_speed_percent=20,
            remove_surface_soak=False,
            min_soak_depth=5,
            max_soak_depth=20,
            use_deck_pressure_offset=False,
            exclude_flags=False,
        )

        expected_flags = expected_data["flag"].values
        assert np.all(expected_flags == result_flags)

    def test_loop_edit_pressure_min_velocity_exclude_flags_pass(self, request):
        expected_data = si.read_cnv_file(
            test_data / "CAST0002_mod_filt_loop_min_v_remove_soak.cnv"
        )
        data = si.read_cnv_file(
            test_data / "CAST0002_mod_filt_loop_min_v_exclude_flags_from_remove_soak.cnv"
        )

        result_flags = sp.loop_edit_pressure(
            pressure=data["prSM"].values,
            latitude=0,
            flag=data["flag"].values,
            sample_interval=data.attrs["sample_interval"],
            min_velocity_type="fixed",
            min_velocity=0.1,
            window_size=3,
            mean_speed_percent=20,
            remove_surface_soak=False,
            min_soak_depth=5,
            max_soak_depth=20,
            use_deck_pressure_offset=False,
            exclude_flags=True,
        )

        expected_flags = expected_data["flag"].values
        mismatches = sum(expected_flags != result_flags)
        assert mismatches == 1

    def test_loop_edit_pressure_mean_speed_percent_remove_soak_pass(self, request):
        expected_data = si.read_cnv_file(
            test_data / "CAST0002_mod_filt_loop_percent_remove_soak.cnv"
        )
        data = si.read_cnv_file(test_data / "CAST0002_mod_filt.cnv")

        result_flags = sp.loop_edit_pressure(
            pressure=data["prSM"].values,
            latitude=0,
            flag=data["flag"].values,
            sample_interval=data.attrs["sample_interval"],
            min_velocity_type="percent",
            min_velocity=0.1,
            window_size=3,
            mean_speed_percent=20,
            remove_surface_soak=True,
            min_soak_depth=5,
            max_soak_depth=20,
            use_deck_pressure_offset=False,
            exclude_flags=False,
        )

        expected_flags = expected_data["flag"].values
        mismatches = sum(expected_flags != result_flags)
        # Less than 2% mismatched flags, outliers appear to conform to spec.
        assert mismatches == 8

    def test_loop_edit_pressure_mean_speed_percent_pass(self, request):
        expected_data = si.read_cnv_file(test_data / "CAST0002_mod_filt_loop_percent.cnv")
        data = si.read_cnv_file(test_data / "CAST0002_mod_filt.cnv")

        result_flags = sp.loop_edit_pressure(
            pressure=data["prSM"].values,
            latitude=0,
            flag=data["flag"].values,
            sample_interval=data.attrs["sample_interval"],
            min_velocity_type="percent",
            min_velocity=0.1,
            window_size=3,
            mean_speed_percent=20,
            remove_surface_soak=False,
            min_soak_depth=5,
            max_soak_depth=20,
            use_deck_pressure_offset=False,
            exclude_flags=False,
        )[100:]

        # SeaSoft is including earliest samples as local maxima during the first
        # downcast which would be discarded when requiring a minimum soak depth
        expected_flags = expected_data["flag"].values[100:]
        mismatches = sum(expected_flags != result_flags)
        # Less than 2% mismatched flags, outliers appear to conform to spec.
        assert mismatches == 6

    def test_loop_edit_pressure_min_velocity_pass_2(self, request):
        expected_data = si.read_cnv_file(test_data / "SBE19plus_loop_edit.cnv")
        data = si.read_cnv_file(test_data / "SBE19plus.cnv")

        result_flags = sp.loop_edit_pressure(
            pressure=data["prdM"].values,
            latitude=0,
            flag=data["flag"].values,
            sample_interval=data.attrs["sample_interval"],
            min_velocity_type="fixed",
            min_velocity=0.25,
            window_size=3,
            mean_speed_percent=20,
            remove_surface_soak=True,
            min_soak_depth=5,
            max_soak_depth=20,
            use_deck_pressure_offset=True,
            exclude_flags=True,
        )[100:]

        # SeaSoft is including earliest samples as local maxima during the first
        # downcast which would be discarded when requiring a minimum soak depth
        expected_flags = expected_data["flag"].values[100:]
        mismatches = sum(expected_flags != result_flags)
        # Less than 2% mismatched flags, outliers appear to conform to spec.
        assert mismatches == 2


class TestBinAverage:
    # fmt: off
    # I think the extra wide lines here are less annoying than the
    # extra tall lines formatted by black. Feel free to enable black
    # here if you disagree

    def test_bin_average_interp_bin_3(self, request):
        source_in = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped.cnv"
        source_out = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped_binavg_interp_bin3.cnv"
        result = si.read_cnv_file(source_in)
        expected = si.read_cnv_file(source_out)
        binavg = sp.bin_average(
            dataset=result,
            bin_variable='prdM',
            bin_size=3,
            interpolate=True
            )
        for variable in list(expected.data_vars):
            tolerance = get_tolerance(expected[variable].values, flag_value=sp.FLAG_VALUE)
            assert np.allclose(expected[variable].values, binavg[variable].values, rtol=0, atol=tolerance)

    def test_bin_average_interp_bin_5(self, request):
        source_in = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped.cnv"
        source_out = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped_binavg_interp_bin5.cnv"
        result = si.read_cnv_file(source_in)
        expected = si.read_cnv_file(source_out)
        binavg = sp.bin_average(
            dataset=result,
            bin_variable='prdM',
            bin_size=5,
            interpolate=True
        )
        for variable in list(expected.data_vars):
            tolerance = get_tolerance(expected[variable].values, flag_value=sp.FLAG_VALUE)
            assert np.allclose(expected[variable].values, binavg[variable].values, rtol=0, atol=tolerance)

    def test_bin_average_default(self, request):
        source_in = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped.cnv"
        source_in2 = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped2.cnv"
        source_out = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped_binavg.cnv"
        source_out2 = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped2_binavg.cnv"
        result = si.read_cnv_file(source_in)
        data2 = si.read_cnv_file(source_in2)
        expected = si.read_cnv_file(source_out)
        expected2 = si.read_cnv_file(source_out2)
        binavg = sp.bin_average(
            dataset = result,
            bin_variable = 'prdM',
            bin_size = 2,
        )
        binavg2 = sp.bin_average(
            dataset = data2,
            bin_variable = 'prdM',
            bin_size = 2,
        )
        for variable in list(expected.data_vars):
            tolerance = get_tolerance(expected[variable].values, flag_value=sp.FLAG_VALUE)
            assert np.allclose(expected[variable].values, binavg[variable].values, rtol=0, atol=tolerance)
            assert np.allclose(expected2[variable], binavg2[variable], rtol=0, atol=tolerance)

    def test_bin_average_include_surface(self, request):
        source_in = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped.cnv"
        source_out = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped_binavg_surface.cnv"
        result = si.read_cnv_file(source_in)
        expected = si.read_cnv_file(source_out)
        binavg = sp.bin_average(
            dataset = result,
            bin_variable = 'prdM',
            bin_size = 2,
            include_surface_bin = True,
            surface_bin_min = 0,
            surface_bin_max = 1,
        )
        for variable in list(expected.data_vars):
            tolerance = get_tolerance(expected[variable].values, flag_value=sp.FLAG_VALUE)
            assert np.allclose(expected[variable].values, binavg[variable].values, rtol=0, atol=tolerance)

    def test_bin_average_interpolate(self, request):
        source_in = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped.cnv"
        source_out = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped_binavg_interp.cnv"
        result = si.read_cnv_file(source_in)
        expected = si.read_cnv_file(source_out)
        binavg = sp.bin_average(
            dataset = result,
            bin_variable = 'prdM',
            bin_size = 2,
            interpolate = True,
        )
        for variable in list(expected.data_vars):
            tolerance = get_tolerance(expected[variable].values, flag_value=sp.FLAG_VALUE)
            assert np.allclose(expected[variable].values, binavg[variable].values, rtol=0, atol=tolerance)

    def test_bin_average_interpolate_surface(self, request):
        source_in = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped.cnv"
        source_out = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped_binavg_interp_surface.cnv"
        result = si.read_cnv_file(source_in)
        expected = si.read_cnv_file(source_out)
        binavg = sp.bin_average(
            dataset = result,
            bin_variable = 'prdM',
            bin_size = 2,
            interpolate = True,
            include_surface_bin = True,
            surface_bin_min = 0,
            surface_bin_max = 1,
            surface_bin_value = 0.5,
        )
        for variable in list(expected.data_vars):
            tolerance = get_tolerance(expected[variable].values, flag_value=sp.FLAG_VALUE)
            assert np.allclose(expected[variable].values, binavg[variable].values, rtol=0, atol=tolerance)

    def test_bin_average_downcast(self, request):
        source_in = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped.cnv"
        source_out = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped_binavg_downcast.cnv"
        result = si.read_cnv_file(source_in)
        expected = si.read_cnv_file(source_out)
        binavg = sp.bin_average(
            dataset = result,
            bin_variable = 'prdM',
            bin_size = 2,
            cast_type = sp.CastType.DOWNCAST
        )
        for variable in list(expected.data_vars):
            tolerance = get_tolerance(expected[variable].values, flag_value=sp.FLAG_VALUE)
            assert np.allclose(expected[variable].values, binavg[variable].values, rtol=0, atol=tolerance)

    def test_bin_average_downcast_interp(self, request):
        source_in = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped.cnv"
        source_out = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped_binavg_downcast_interp.cnv"
        result = si.read_cnv_file(source_in)
        expected = si.read_cnv_file(source_out)
        binavg = sp.bin_average(
            dataset = result,
            bin_variable = 'prdM',
            bin_size = 2,
            cast_type = sp.CastType.DOWNCAST,
            interpolate = True,
        )
        for variable in list(expected.data_vars):
            tolerance = get_tolerance(expected[variable].values, flag_value=sp.FLAG_VALUE)
            assert np.allclose(expected[variable].values, binavg[variable].values, rtol=0, atol=tolerance)

    def test_bin_average_downcast_interp_surface(self, request):
        source_in = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped.cnv"
        source_out = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped_binavg_downcast_interp_surface.cnv"
        result = si.read_cnv_file(source_in)
        expected = si.read_cnv_file(source_out)
        binavg = sp.bin_average(
            dataset = result,
            bin_variable = 'prdM',
            bin_size = 2,
            cast_type = sp.CastType.DOWNCAST,
            interpolate = True,
            include_surface_bin = True,
            surface_bin_min = 0,
            surface_bin_max = 1,
            surface_bin_value = 0.5,
        )
        for variable in list(expected.data_vars):
            tolerance = get_tolerance(expected[variable].values, flag_value=sp.FLAG_VALUE)
            assert np.allclose(expected[variable].values, binavg[variable].values, rtol=0, atol=tolerance)

    def test_bin_average_upcast(self, request):
        source_in = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped.cnv"
        source_out = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped_binavg_upcast.cnv"
        result = si.read_cnv_file(source_in)
        expected = si.read_cnv_file(source_out)
        binavg = sp.bin_average(
            dataset = result,
            bin_variable = 'prdM',
            bin_size = 2,
            cast_type = sp.CastType.UPCAST
        )
        for variable in list(expected.data_vars):
            tolerance = get_tolerance(expected[variable].values, flag_value=sp.FLAG_VALUE)
            assert np.allclose(expected[variable].values, binavg[variable].values, rtol=0, atol=tolerance)

    def test_bin_average_upcast_interp(self, request):
        source_in = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped.cnv"
        source_out = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped_binavg_upcast_interp.cnv"
        result = si.read_cnv_file(source_in)
        expected = si.read_cnv_file(source_out)
        binavg = sp.bin_average(
            dataset = result,
            bin_variable = 'prdM',
            bin_size = 2,
            cast_type = sp.CastType.UPCAST,
            interpolate = True,
        )
        for variable in list(expected.data_vars):
            tolerance = get_tolerance(expected[variable].values, flag_value=sp.FLAG_VALUE)
            assert np.allclose(expected[variable].values, binavg[variable].values, rtol=0, atol=tolerance)

    def test_bin_average_upcast_interp_surface(self, request):
        source_in = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped.cnv"
        source_out = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped_binavg_upcast_interp_surface.cnv"
        result = si.read_cnv_file(source_in)
        expected = si.read_cnv_file(source_out)
        binavg = sp.bin_average(
            dataset = result,
            bin_variable = 'prdM',
            bin_size = 2,
            cast_type = sp.CastType.UPCAST,
            interpolate = True,
            include_surface_bin = True,
            surface_bin_min = 0,
            surface_bin_max = 1,
            surface_bin_value = 0.5,
        )
        for variable in list(expected.data_vars):
            tolerance = get_tolerance(expected[variable].values, flag_value=sp.FLAG_VALUE)
            assert np.allclose(expected[variable].values, binavg[variable].values, rtol=0, atol=tolerance)

    def test_bin_average_exclude_flags(self, request):
        source_in = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped_loopedit_wildedit.cnv"
        source_out = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped_loopedit_wildedit_binavg_exclude.cnv"
        result = si.read_cnv_file(source_in)
        expected = si.read_cnv_file(source_out)
        binavg = sp.bin_average(
            dataset = result,
            bin_variable = 'prdM',
            bin_size = 2,
        )
        for variable in list(expected.data_vars):
            tolerance = get_tolerance(expected[variable].values, flag_value=sp.FLAG_VALUE)
            assert np.allclose(expected[variable].values, binavg[variable].values, rtol=0, atol=tolerance)

    def test_bin_average_include_flags(self, request):
        source_in = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped_loopedit_wildedit.cnv"
        source_out = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped_loopedit_wildedit_binavg_include.cnv"
        result = si.read_cnv_file(source_in)
        expected = si.read_cnv_file(source_out)
        binavg = sp.bin_average(
            dataset = result,
            bin_variable = 'prdM',
            bin_size = 2,
            exclude_bad_scans = False
        )
        for variable in list(expected.data_vars):
            tolerance = get_tolerance(expected[variable].values, flag_value=sp.FLAG_VALUE)
            assert np.allclose(expected[variable].values, binavg[variable].values, rtol=0, atol=tolerance)

    def test_bin_average_include_flags_interp(self, request):
        source_in = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped_loopedit_wildedit.cnv"
        source_out = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped_loopedit_wildedit_binavg_include_interp.cnv"
        result = si.read_cnv_file(source_in)
        expected = si.read_cnv_file(source_out)
        binavg = sp.bin_average(
            dataset = result,
            bin_variable = 'prdM',
            bin_size = 2,
            interpolate = True,
            exclude_bad_scans = False,
        )
        for variable in list(expected.data_vars):
            tolerance = get_tolerance(expected[variable].values, flag_value=sp.FLAG_VALUE)
            assert np.allclose(expected[variable].values, binavg[variable].values, rtol=0, atol=tolerance)

    def test_bin_average_exclude_flags_interp(self, request):
        source_in = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped_loopedit_wildedit.cnv"
        source_out = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped_loopedit_wildedit_binavg_exclude_interp.cnv"
        result = si.read_cnv_file(source_in)
        expected = si.read_cnv_file(source_out)
        binavg = sp.bin_average(
            dataset = result,
            bin_variable = 'prdM',
            bin_size = 2,
            interpolate = True,
        )
        for variable in list(expected.data_vars):
            tolerance = get_tolerance(expected[variable].values, flag_value=sp.FLAG_VALUE)
            assert np.allclose(expected[variable].values, binavg[variable].values, rtol=0, atol=tolerance)

    def test_bin_average_time(self, request):
        source_in = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped.cnv"
        source_out = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped_binavg_time.cnv"
        result = si.read_cnv_file(source_in)
        expected = si.read_cnv_file(source_out)
        binavg = sp.bin_average(
            dataset = result,
            bin_variable = 'timeS',
            bin_size = 10,
            cast_type = sp.CastType.NONE
        )
        for variable in list(expected.data_vars):
            tolerance = get_tolerance(expected[variable].values, flag_value=sp.FLAG_VALUE)
            assert np.allclose(expected[variable].values, binavg[variable].values, rtol=0, atol=tolerance)

    def test_bin_average_scan(self, request):
        source_in = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped.cnv"
        source_out = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped_binavg_scan.cnv"
        result = si.read_cnv_file(source_in)
        expected = si.read_cnv_file(source_out)
        binavg = sp.bin_average(
            dataset = result,
            bin_variable = 'nScan',
            bin_size = 10,
            cast_type = sp.CastType.NONE
        )
        for variable in list(expected.data_vars):
            tolerance = get_tolerance(expected[variable].values, flag_value=sp.FLAG_VALUE)
            assert np.allclose(expected[variable].values, binavg[variable].values, rtol=0, atol=tolerance)

    # fmt: on


class TestWildEdit:
    def test_wild_edit_pass(self, request):
        expected_dataset = si.read_cnv_file(
            test_data / "19plus_V2_CTD-processing_example_wild_edit.cnv"
        )
        expected_conductivity = expected_dataset["c0S/m"].values

        dataset = si.read_cnv_file(test_data / "19plus_V2_CTD-processing_example.cnv")
        conductivity = dataset["c0S/m"].values
        flags = dataset["flag"].values

        wild_edit_output = sp.wild_edit(conductivity, flags, 2, 20, 100, 0, False)

        request.node.return_value = wild_edit_output.tolist()
        assert np.all(wild_edit_output == expected_conductivity)


class TestWindowFilter:
    file_prefix = Path(test_data / "19plus_V2_CTD-processing_example")
    # cnvdata = si.read_cnv_file(f"{file_prefix}.cnv")
    # pressure = cnvdata["prdM"].values
    # flags = cnvdata["flag"].values
    window_width = 5
    half_width = 1  # only applies to gaussian
    offset = 0.25  # only applies to gaussian
    # sample_interval = cnvdata.attrs["sample_interval"]  # only applies to gaussian

    @pytest.fixture
    def cnvdata(self):
        return si.read_cnv_file(f"{self.file_prefix}.cnv")

    def test_boxcar_filter(self, cnvdata, request):
        expected_dataset = si.read_cnv_file(f"{self.file_prefix}_boxcar_5.cnv")
        expected_pressure = expected_dataset["prdM"].values

        filtered_pressure = sp.window_filter(
            cnvdata["prdM"].values,
            cnvdata["flag"].values,
            "boxcar",
            self.window_width,
            cnvdata.attrs["sample_interval"],
            self.half_width,
            self.offset,
            False,
        )

        request.node.return_value = filtered_pressure.tolist()
        assert close_enough(filtered_pressure, expected_pressure, 3, 1e-12)

    def test_boxcar_filter_exclude_flags(self, cnvdata, request):
        expected_dataset = si.read_cnv_file(f"{self.file_prefix}_boxcar_5_excluded.cnv")
        expected_pressure = expected_dataset["prdM"].values

        filtered_pressure = sp.window_filter(
            cnvdata["prdM"].values,
            cnvdata["flag"].values,
            "boxcar",
            self.window_width,
            cnvdata.attrs["sample_interval"],
            self.half_width,
            self.offset,
            True,
        )

        request.node.return_value = filtered_pressure.tolist()
        assert close_enough(filtered_pressure, expected_pressure, 3, 1e-12)

    def test_cosine_filter(self, cnvdata, request):
        expected_dataset = si.read_cnv_file(f"{self.file_prefix}_cosine_5.cnv")
        expected_pressure = expected_dataset["prdM"].values

        filtered_pressure = sp.window_filter(
            cnvdata["prdM"].values,
            cnvdata["flag"].values,
            "cosine",
            self.window_width,
            cnvdata.attrs["sample_interval"],
            self.half_width,
            self.offset,
            False,
        )

        request.node.return_value = filtered_pressure.tolist()
        assert close_enough(filtered_pressure, expected_pressure, 3, 1e-12)

    def test_cosine_filter_exclude_flags(self, cnvdata, request):
        expected_dataset = si.read_cnv_file(f"{self.file_prefix}_cosine_5_excluded.cnv")
        expected_pressure = expected_dataset["prdM"].values

        filtered_pressure = sp.window_filter(
            cnvdata["prdM"].values,
            cnvdata["flag"].values,
            "cosine",
            self.window_width,
            cnvdata.attrs["sample_interval"],
            self.half_width,
            self.offset,
            True,
        )

        request.node.return_value = filtered_pressure.tolist()
        assert close_enough(filtered_pressure, expected_pressure, 3, 1e-12)

    def test_triangle_filter(self, cnvdata, request):
        expected_dataset = si.read_cnv_file(f"{self.file_prefix}_triangle_5.cnv")
        expected_pressure = expected_dataset["prdM"].values

        filtered_pressure = sp.window_filter(
            cnvdata["prdM"].values,
            cnvdata["flag"].values,
            "triangle",
            self.window_width,
            cnvdata.attrs["sample_interval"],
            self.half_width,
            self.offset,
            False,
        )

        request.node.return_value = filtered_pressure.tolist()
        assert close_enough(filtered_pressure, expected_pressure, 3, 1e-12)

    def test_triangle_filter_exclude_flags(self, cnvdata, request):
        expected_dataset = si.read_cnv_file(f"{self.file_prefix}_triangle_5_excluded.cnv")
        expected_pressure = expected_dataset["prdM"].values

        filtered_pressure = sp.window_filter(
            cnvdata["prdM"].values,
            cnvdata["flag"].values,
            "triangle",
            self.window_width,
            cnvdata.attrs["sample_interval"],
            self.half_width,
            self.offset,
            True,
        )

        request.node.return_value = filtered_pressure.tolist()
        assert close_enough(filtered_pressure, expected_pressure, 3, 1e-12)

    def test_gaussian_filter(self, cnvdata, request):
        expected_dataset = si.read_cnv_file(f"{self.file_prefix}_gaussian_5_1_025.cnv")
        expected_pressure = expected_dataset["prdM"].values

        filtered_pressure = sp.window_filter(
            cnvdata["prdM"].values,
            cnvdata["flag"].values,
            "gaussian",
            self.window_width,
            cnvdata.attrs["sample_interval"],
            self.half_width,
            self.offset,
            False,
        )

        request.node.return_value = filtered_pressure.tolist()
        assert close_enough(filtered_pressure, expected_pressure, 3, 1e-12)

    def test_gaussian_filter_exclude_flags(self, cnvdata, request):
        expected_dataset = si.read_cnv_file(f"{self.file_prefix}_gaussian_5_1_025_excluded.cnv")
        expected_pressure = expected_dataset["prdM"].values

        filtered_pressure = sp.window_filter(
            cnvdata["prdM"].values,
            cnvdata["flag"].values,
            "gaussian",
            self.window_width,
            cnvdata.attrs["sample_interval"],
            self.half_width,
            self.offset,
            True,
        )

        request.node.return_value = filtered_pressure.tolist()
        assert close_enough(filtered_pressure, expected_pressure, 3, 1e-12)

    def test_median_filter(self, cnvdata, request):
        expected_dataset = si.read_cnv_file(f"{self.file_prefix}_median_5.cnv")
        expected_pressure = expected_dataset["prdM"].values

        filtered_pressure = sp.window_filter(
            cnvdata["prdM"].values,
            cnvdata["flag"].values,
            "median",
            self.window_width,
            cnvdata.attrs["sample_interval"],
            self.half_width,
            self.offset,
            False,
        )
        request.node.return_value = filtered_pressure.tolist()
        assert close_enough(filtered_pressure, expected_pressure, 3, 1e-12)

    def test_median_filter_exclude_flags(self, cnvdata, request):
        expected_dataset = si.read_cnv_file(f"{self.file_prefix}_median_5_excluded.cnv")
        expected_pressure = expected_dataset["prdM"].values

        filtered_pressure = sp.window_filter(
            cnvdata["prdM"].values,
            cnvdata["flag"].values,
            "median",
            self.window_width,
            cnvdata.attrs["sample_interval"],
            self.half_width,
            self.offset,
            True,
        )
        request.node.return_value = filtered_pressure.tolist()
        assert close_enough(filtered_pressure, expected_pressure, 3, 1e-12)


class TestNitrate:
    dac_min = -5
    dac_max = 100
    voltages = np.array([0.095, 1, 2, 3, 4.095])

    def test_convert_nitrate_umno3(self, request):
        expected_umno3 = np.array([-5, 18.75625, 45.00625, 71.25625, 100])
        nitrate = sc.convert_nitrate(self.voltages, dac_min=self.dac_min, dac_max=self.dac_max)
        request.node.return_value = nitrate.tolist()
        assert np.allclose(expected_umno3, nitrate, atol=0.000001)

    def test_convert_nitrate_mgnl(self, request):
        expected_mgnl = np.array([-0.070035, 0.26271879375, 0.63040254375, 0.99808629375, 1.4007])
        nitrate = sc.convert_nitrate(
            self.voltages, dac_min=self.dac_min, dac_max=self.dac_max, units="mgNL"
        )
        request.node.return_value = nitrate.tolist()
        assert np.allclose(expected_mgnl, nitrate, atol=0.000001)


class TestSplit:
    @pytest.mark.parametrize(
        "source_path, cast_type",
        [
            ("SBE19plus_01906398_2019_07_15_0033_cropped_upcast.cnv", sp.CastType.UPCAST),
            ("SBE19plus_01906398_2019_07_15_0033_cropped_downcast.cnv", sp.CastType.DOWNCAST),
        ],
    )
    def test_split(self, source_path, cast_type, request):
        expected_source = test_data / source_path
        expected = si.read_cnv_file(expected_source)

        source = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped.cnv"
        result = si.read_cnv_file(source)
        result = sp.split(result, "prdM", cast_type=cast_type, drop=True)

        for variable in list(expected.data_vars):
            tolerance = get_tolerance(expected[variable].values, flag_value=sp.FLAG_VALUE)
            assert np.allclose(
                expected[variable].values, result[variable].values, rtol=0, atol=10**tolerance
            )

    @pytest.mark.parametrize(
        "source_path, cast_type",
        [
            ("SBE19plus_01906398_2019_07_15_0033_cropped_loopedit_upcast.cnv", sp.CastType.UPCAST),
            (
                "SBE19plus_01906398_2019_07_15_0033_cropped_loopedit_downcast.cnv",
                sp.CastType.DOWNCAST,
            ),
        ],
    )
    def test_split_with_flags(self, source_path, cast_type, request):
        expected_source = test_data / source_path
        expected = si.read_cnv_file(expected_source)

        source = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped_loopedit.cnv"
        result = si.read_cnv_file(source)
        result = sp.split(result, "prdM", exclude_bad_scans=True, cast_type=cast_type, drop=True)

        for variable in list(expected.data_vars):
            tolerance = get_tolerance(expected[variable].values, flag_value=sp.FLAG_VALUE)
            assert np.allclose(
                expected[variable].values, result[variable].values, rtol=0, atol=10**tolerance
            )
