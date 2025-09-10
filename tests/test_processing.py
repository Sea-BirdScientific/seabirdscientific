"""Data processing unit tests.

"""

# Native imports
from logging import getLogger
from pathlib import Path

# Third-party imports
import numpy as np
import pandas as pd
import pytest

# Sea-Bird imports

# Internal imports
import seabirdscientific.instrument_data as idata
import seabirdscientific.processing as p
import seabirdscientific.conversion as c
from seabirdscientific.utils import close_enough, get_decimal_length

test_data = Path("./tests/resources/test-data")
logger = getLogger(__name__)


class TestLowPassFilter:

    def test_low_pass_filter(self, request):
        expected_path = test_data / "SBE37SM-filtered.asc"
        expected = pd.read_csv(expected_path)["Tv290C"].values

        source_path = test_data / "SBE37SM-unfiltered.asc"
        source = pd.read_csv(source_path)["Tv290C"].values

        filtered = p.low_pass_filter(source, 10000, sample_interval=120)

        request.node.return_value = filtered.tolist()

        assert np.allclose(expected, filtered, rtol=0, atol=1e-4)


class TestAlignCtd:
    expected_data_path = test_data / "SBE37SM-align.cnv"
    expected_data = idata.cnv_to_instrument_data(expected_data_path)
    # expected data was algined with:
    #   sample_interval=120,
    #   tv290C + 12s
    #   cond0S/m + 150s
    #   prdM + 120s
    #   prdE + 0s
    #   sal00 - 240s

    source_data_path = test_data / "SBE37SM.cnv"
    source_data = idata.cnv_to_instrument_data(source_data_path)

    def test_align_ctd_add_simple_pass(self, request):
        # Fix last valid sample
        self.expected_data.measurements["tv290C"].values[
            len(self.expected_data.measurements["tv290C"].values) - 1
        ] = -9.99e-29
        result = p.align_ctd(self.source_data.measurements["tv290C"].values, 12, 120)
        request.node.return_value = result.tolist()
        assert np.allclose(self.expected_data.measurements["tv290C"].values, result, atol=0.0001)

    def test_align_ctd_add_simple_pass_2(self, request):
        # Fix last valid sample
        self.expected_data.measurements["cond0S/m"].values[
            len(self.expected_data.measurements["cond0S/m"].values) - 2
        ] = -9.99e-29
        result = p.align_ctd(self.source_data.measurements["cond0S/m"].values, 150, 120)
        request.node.return_value = result.tolist()
        assert np.allclose(self.expected_data.measurements["cond0S/m"].values, result, atol=0.0001)

    def test_align_ctd_add_exact_factor(self, request):
        result = p.align_ctd(self.source_data.measurements["prdM"].values, 120, 120)
        request.node.return_value = result.tolist()
        assert np.allclose(self.expected_data.measurements["prdM"].values, result, atol=0.0001)

    def test_align_ctd_no_change(self, request):
        result = p.align_ctd(self.source_data.measurements["prdE"].values, 0, 120)
        request.node.return_value = result.tolist()
        assert np.allclose(self.expected_data.measurements["prdE"].values, result, atol=0.0001)

    def test_align_ctd_subtract(self, request):
        result = p.align_ctd(self.source_data.measurements["sal00"].values, -240, 120)
        request.node.return_value = result.tolist()
        assert np.allclose(self.expected_data.measurements["sal00"].values, result, atol=0.0001)


class TestCellThermalMass:
    def test_cell_thermal_mass_pass(self, request):
        expected_data_path = test_data / "SBE37SM-ctm.cnv"
        expected_data = idata.cnv_to_instrument_data(expected_data_path)
        source_data_path = test_data / "SBE37SM.cnv"
        source_data = idata.cnv_to_instrument_data(source_data_path)
        corrected_conductivity = p.cell_thermal_mass(
            source_data.measurements["tv290C"].values,
            source_data.measurements["cond0S/m"].values,
            0.03,
            7.0,
            120.0,
        )
        request.node.return_value = corrected_conductivity.tolist()
        assert np.allclose(
            expected_data.measurements["cond0S/m"].values,
            corrected_conductivity,
            atol=0.000001,
        )


class TestLoopEdit:
    def test_loop_edit_pressure_min_velocity_pass(self, request):
        expected_data = idata.cnv_to_instrument_data(
            test_data / "CAST0002_mod_filt_loop_min_v.cnv"
        )
        data = idata.cnv_to_instrument_data(test_data / "CAST0002_mod_filt.cnv")

        p.loop_edit_pressure(
            pressure=data.measurements["prSM"].values,
            latitude=data.latitude,
            flag=data.measurements["flag"].values,
            sample_interval=data.interval_s,
            min_velocity_type=p.MinVelocityType.FIXED,
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

        request.node.return_value = result_flags.tolist()
        assert mismatches / len(result_flags) < 0.01

    def test_loop_edit_pressure_min_velocity_remove_soak_pass(self, request):
        expected_data = idata.cnv_to_instrument_data(
            test_data / "CAST0002_mod_filt_loop_min_v_remove_soak.cnv"
        )
        data = idata.cnv_to_instrument_data(test_data / "CAST0002_mod_filt.cnv")

        p.loop_edit_pressure(
            pressure=data.measurements["prSM"].values,
            latitude=data.latitude,
            flag=data.measurements["flag"].values,
            sample_interval=data.interval_s,
            min_velocity_type=p.MinVelocityType.FIXED,
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

        request.node.return_value = result_flags.tolist()
        assert mismatches / len(result_flags) < 0.01

    def test_loop_edit_pressure_min_velocity_reset_flags_pass(self, request):
        expected_data = idata.cnv_to_instrument_data(
            test_data / "CAST0002_mod_filt_loop_min_v.cnv"
        )
        data = idata.cnv_to_instrument_data(
            test_data / "CAST0002_mod_filt_loop_min_v_remove_soak.cnv"
        )

        p.loop_edit_pressure(
            pressure=data.measurements["prSM"].values,
            latitude=data.latitude,
            flag=data.measurements["flag"].values,
            sample_interval=data.interval_s,
            min_velocity_type=p.MinVelocityType.FIXED,
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

        request.node.return_value = result_flags.tolist()
        assert mismatches / len(result_flags) < 0.01

    def test_loop_edit_pressure_min_velocity_exclude_flags_pass(self, request):
        expected_data = idata.cnv_to_instrument_data(
            test_data / "CAST0002_mod_filt_loop_min_v_remove_soak.cnv"
        )
        data = idata.cnv_to_instrument_data(
            test_data / "CAST0002_mod_filt_loop_min_v_exclude_flags_from_remove_soak.cnv"
        )

        p.loop_edit_pressure(
            pressure=data.measurements["prSM"].values,
            latitude=data.latitude,
            flag=data.measurements["flag"].values,
            sample_interval=data.interval_s,
            min_velocity_type=p.MinVelocityType.FIXED,
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

        request.node.return_value = result_flags.tolist()
        assert mismatches / len(result_flags) < 0.01

    def test_loop_edit_pressure_mean_speed_percent_remove_soak_pass(self, request):
        expected_data = idata.cnv_to_instrument_data(
            test_data / "CAST0002_mod_filt_loop_percent_remove_soak.cnv"
        )
        data = idata.cnv_to_instrument_data(test_data / "CAST0002_mod_filt.cnv")

        p.loop_edit_pressure(
            pressure=data.measurements["prSM"].values,
            latitude=data.latitude,
            flag=data.measurements["flag"].values,
            sample_interval=data.interval_s,
            min_velocity_type=p.MinVelocityType.PERCENT,
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

        request.node.return_value = result_flags.tolist()
        # Less than 2% mismatched flags, outliers appear to conform to spec.
        assert mismatches / len(result_flags) < 0.02

    def test_loop_edit_pressure_mean_speed_percent_pass(self, request):
        expected_data = idata.cnv_to_instrument_data(
            test_data / "CAST0002_mod_filt_loop_percent.cnv"
        )
        data = idata.cnv_to_instrument_data(test_data / "CAST0002_mod_filt.cnv")

        p.loop_edit_pressure(
            pressure=data.measurements["prSM"].values,
            latitude=data.latitude,
            flag=data.measurements["flag"].values,
            sample_interval=data.interval_s,
            min_velocity_type=p.MinVelocityType.PERCENT,
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

        request.node.return_value = result_flags.tolist()
        # Less than 2% mismatched flags, outliers appear to conform to spec.
        assert mismatches / len(result_flags) < 0.02

    def test_loop_edit_pressure_min_velocity_pass_2(self, request):
        expected_data = idata.cnv_to_instrument_data(
            test_data / "SBE19plus_loop_edit.cnv"
        )
        data = idata.cnv_to_instrument_data(test_data / "SBE19plus.cnv")

        p.loop_edit_pressure(
            pressure=data.measurements["prdM"].values,
            latitude=data.latitude,
            flag=data.measurements["flag"].values,
            sample_interval=data.interval_s,
            min_velocity_type=p.MinVelocityType.FIXED,
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

        request.node.return_value = result_flags.tolist()
        # Less than 2% mismatched flags, outliers appear to conform to spec.
        assert mismatches / len(result_flags) < 0.02


class TestBinAverage:
    # fmt: off
    # I think the extra wide lines here are less annoying than the
    # extra tall lines formatted by black. Feel free to enable black
    # here if you disagree

    def test_bin_average_interp_bin_3(self, request):
        source_in = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped.cnv"
        source_out = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped_binavg_interp_bin3.cnv"
        data = idata.cnv_to_instrument_data(source_in)._to_dataframe()
        expected = idata.cnv_to_instrument_data(source_out)._to_dataframe()
        binavg = p.bin_average(
            dataset=data,
            bin_variable='prdM',
            bin_size=3,
            interpolate=True
            )
        for variable in expected.columns:
            exponent = -1 * get_decimal_length(expected[variable])
            assert np.allclose(expected[variable], binavg[variable], rtol=0, atol=10**exponent)

    def test_bin_average_interp_bin_5(self, request):
        source_in = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped.cnv"
        source_out = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped_binavg_interp_bin5.cnv"
        data = idata.cnv_to_instrument_data(source_in)._to_dataframe()
        expected = idata.cnv_to_instrument_data(source_out)._to_dataframe()
        binavg = p.bin_average(
            dataset=data,
            bin_variable='prdM',
            bin_size=5,
            interpolate=True
        )
        for variable in expected.columns:
            exponent = -1 * get_decimal_length(expected[variable])
            assert np.allclose(expected[variable], binavg[variable], rtol=0, atol=10**exponent)

    def test_bin_average_default(self, request):
        source_in = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped.cnv"
        source_in2 = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped2.cnv"
        source_out = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped_binavg.cnv"
        source_out2 = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped2_binavg.cnv"
        data = idata.cnv_to_instrument_data(source_in)._to_dataframe()
        data2 = idata.cnv_to_instrument_data(source_in2)._to_dataframe()
        expected = idata.cnv_to_instrument_data(source_out)._to_dataframe()
        expected2 = idata.cnv_to_instrument_data(source_out2)._to_dataframe()
        binavg = p.bin_average(
            dataset = data,
            bin_variable = 'prdM',
            bin_size = 2,
        )
        binavg2 = p.bin_average(
            dataset = data2,
            bin_variable = 'prdM',
            bin_size = 2,
        )
        for variable in expected.columns:
            exponent = -1 * get_decimal_length(expected[variable])
            assert np.allclose(expected[variable], binavg[variable], rtol=0, atol=10**exponent)
            assert np.allclose(expected2[variable], binavg2[variable], rtol=0, atol=10**exponent)

    def test_bin_average_include_surface(self, request):
        source_in = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped.cnv"
        source_out = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped_binavg_surface.cnv"
        data = idata.cnv_to_instrument_data(source_in)._to_dataframe()
        expected = idata.cnv_to_instrument_data(source_out)._to_dataframe()
        binavg = p.bin_average(
            dataset = data,
            bin_variable = 'prdM',
            bin_size = 2,
            include_surface_bin = True,
            surface_bin_min = 0,
            surface_bin_max = 1,
        )
        for variable in expected.columns:
            exponent = -1 * get_decimal_length(expected[variable])
            assert np.allclose(expected[variable], binavg[variable], rtol=0, atol=10**exponent)

    def test_bin_average_interpolate(self, request):
        source_in = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped.cnv"
        source_out = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped_binavg_interp.cnv"
        data = idata.cnv_to_instrument_data(source_in)._to_dataframe()
        expected = idata.cnv_to_instrument_data(source_out)._to_dataframe()
        binavg = p.bin_average(
            dataset = data,
            bin_variable = 'prdM',
            bin_size = 2,
            interpolate = True,
        )
        for variable in expected.columns:
            exponent = -1 * get_decimal_length(expected[variable])
            assert np.allclose(expected[variable], binavg[variable], rtol=0, atol=10**exponent)
    
    def test_bin_average_interpolate_surface(self, request):
        source_in = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped.cnv"
        source_out = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped_binavg_interp_surface.cnv"
        data = idata.cnv_to_instrument_data(source_in)._to_dataframe()
        expected = idata.cnv_to_instrument_data(source_out)._to_dataframe()
        binavg = p.bin_average(
            dataset = data,
            bin_variable = 'prdM',
            bin_size = 2,
            interpolate = True,
            include_surface_bin = True,
            surface_bin_min = 0,
            surface_bin_max = 1,
            surface_bin_value = 0.5,
        )
        for variable in expected.columns:
            exponent = -1 * get_decimal_length(expected[variable])
            assert np.allclose(expected[variable], binavg[variable], rtol=0, atol=10**exponent)

    def test_bin_average_downcast(self, request):
        source_in = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped.cnv"
        source_out = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped_binavg_downcast.cnv"
        data = idata.cnv_to_instrument_data(source_in)._to_dataframe()
        expected = idata.cnv_to_instrument_data(source_out)._to_dataframe()
        binavg = p.bin_average(
            dataset = data,
            bin_variable = 'prdM',
            bin_size = 2,
            cast_type = p.CastType.DOWNCAST
        )
        for variable in expected.columns:
            exponent = -1 * get_decimal_length(expected[variable])
            assert np.allclose(expected[variable], binavg[variable], rtol=0, atol=10**exponent)

    def test_bin_average_downcast_interp(self, request):
        source_in = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped.cnv"
        source_out = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped_binavg_downcast_interp.cnv"
        data = idata.cnv_to_instrument_data(source_in)._to_dataframe()
        expected = idata.cnv_to_instrument_data(source_out)._to_dataframe()
        binavg = p.bin_average(
            dataset = data,
            bin_variable = 'prdM',
            bin_size = 2,
            cast_type = p.CastType.DOWNCAST,
            interpolate = True,
        )
        for variable in expected.columns:
            exponent = -1 * get_decimal_length(expected[variable])
            assert np.allclose(expected[variable], binavg[variable], rtol=0, atol=10**exponent)

    def test_bin_average_downcast_interp_surface(self, request):
        source_in = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped.cnv"
        source_out = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped_binavg_downcast_interp_surface.cnv"
        data = idata.cnv_to_instrument_data(source_in)._to_dataframe()
        expected = idata.cnv_to_instrument_data(source_out)._to_dataframe()
        binavg = p.bin_average(
            dataset = data,
            bin_variable = 'prdM',
            bin_size = 2,
            cast_type = p.CastType.DOWNCAST,
            interpolate = True,
            include_surface_bin = True,
            surface_bin_min = 0,
            surface_bin_max = 1,
            surface_bin_value = 0.5,
        )
        for variable in expected.columns:
            exponent = -1 * get_decimal_length(expected[variable])
            assert np.allclose(expected[variable], binavg[variable], rtol=0, atol=10**exponent)

    def test_bin_average_upcast(self, request):
        source_in = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped.cnv"
        source_out = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped_binavg_upcast.cnv"
        data = idata.cnv_to_instrument_data(source_in)._to_dataframe()
        expected = idata.cnv_to_instrument_data(source_out)._to_dataframe()
        binavg = p.bin_average(
            dataset = data,
            bin_variable = 'prdM',
            bin_size = 2,
            cast_type = p.CastType.UPCAST
        )
        for variable in expected.columns:
            exponent = -1 * get_decimal_length(expected[variable])
            assert np.allclose(expected[variable], binavg[variable], rtol=0, atol=10**exponent)

    def test_bin_average_upcast_interp(self, request):
        source_in = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped.cnv"
        source_out = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped_binavg_upcast_interp.cnv"
        data = idata.cnv_to_instrument_data(source_in)._to_dataframe()
        expected = idata.cnv_to_instrument_data(source_out)._to_dataframe()
        binavg = p.bin_average(
            dataset = data,
            bin_variable = 'prdM',
            bin_size = 2,
            cast_type = p.CastType.UPCAST,
            interpolate = True,
        )
        for variable in expected.columns:
            exponent = -1 * get_decimal_length(expected[variable])
            assert np.allclose(expected[variable], binavg[variable], rtol=0, atol=10**exponent)

    def test_bin_average_upcast_interp_surface(self, request):
        source_in = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped.cnv"
        source_out = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped_binavg_upcast_interp_surface.cnv"
        data = idata.cnv_to_instrument_data(source_in)._to_dataframe()
        expected = idata.cnv_to_instrument_data(source_out)._to_dataframe()
        binavg = p.bin_average(
            dataset = data,
            bin_variable = 'prdM',
            bin_size = 2,
            cast_type = p.CastType.UPCAST,
            interpolate = True,
            include_surface_bin = True,
            surface_bin_min = 0,
            surface_bin_max = 1,
            surface_bin_value = 0.5,
        )
        for variable in expected.columns:
            exponent = -1 * get_decimal_length(expected[variable])
            assert np.allclose(expected[variable], binavg[variable], rtol=0, atol=10**exponent)

    def test_bin_average_exclude_flags(self, request):
        source_in = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped_loopedit_wildedit.cnv"
        source_out = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped_loopedit_wildedit_binavg_exclude.cnv"
        data = idata.cnv_to_instrument_data(source_in)._to_dataframe()
        expected = idata.cnv_to_instrument_data(source_out)._to_dataframe()
        binavg = p.bin_average(
            dataset = data,
            bin_variable = 'prdM',
            bin_size = 2,
        )
        for variable in expected.columns:
            exponent = -1 * get_decimal_length(expected[variable])
            assert np.allclose(expected[variable], binavg[variable], rtol=0, atol=10**exponent)

    def test_bin_average_include_flags(self, request):
        source_in = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped_loopedit_wildedit.cnv"
        source_out = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped_loopedit_wildedit_binavg_include.cnv"
        data = idata.cnv_to_instrument_data(source_in)._to_dataframe()
        expected = idata.cnv_to_instrument_data(source_out)._to_dataframe()
        binavg = p.bin_average(
            dataset = data,
            bin_variable = 'prdM',
            bin_size = 2,
            exclude_bad_scans = False
        )
        for variable in expected.columns:
            exponent = -1 * get_decimal_length(expected[variable])
            assert np.allclose(expected[variable], binavg[variable], rtol=0, atol=10**exponent)

    def test_bin_average_include_flags_interp(self, request):
        source_in = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped_loopedit_wildedit.cnv"
        source_out = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped_loopedit_wildedit_binavg_include_interp.cnv"
        data = idata.cnv_to_instrument_data(source_in)._to_dataframe()
        expected = idata.cnv_to_instrument_data(source_out)._to_dataframe()
        binavg = p.bin_average(
            dataset = data,
            bin_variable = 'prdM',
            bin_size = 2,
            interpolate = True,
            exclude_bad_scans = False,
        )
        for variable in expected.columns:
            exponent = -1 * get_decimal_length(expected[variable])
            assert np.allclose(expected[variable], binavg[variable], rtol=0, atol=10**exponent)

    def test_bin_average_exclude_flags_interp(self, request):
        source_in = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped_loopedit_wildedit.cnv"
        source_out = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped_loopedit_wildedit_binavg_exclude_interp.cnv"
        data = idata.cnv_to_instrument_data(source_in)._to_dataframe()
        expected = idata.cnv_to_instrument_data(source_out)._to_dataframe()
        binavg = p.bin_average(
            dataset = data,
            bin_variable = 'prdM',
            bin_size = 2,
            interpolate = True,
        )
        for variable in expected.columns:
            exponent = -1 * get_decimal_length(expected[variable])
            assert np.allclose(expected[variable], binavg[variable], rtol=0, atol=10**exponent)

    def test_bin_average_time(self, request):
        source_in = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped.cnv"
        source_out = test_data / "SBE19plus_01906398_2019_07_15_0033_cropped_binavg_time.cnv"
        data = idata.cnv_to_instrument_data(source_in)._to_dataframe()
        expected = idata.cnv_to_instrument_data(source_out)._to_dataframe()
        binavg = p.bin_average(
            dataset = data,
            bin_variable = 'timeS',
            bin_size = 10,
            cast_type = p.CastType.NA
        )
        for variable in expected.columns:
            exponent = -1 * get_decimal_length(expected[variable])
            assert np.allclose(expected[variable], binavg[variable], rtol=0, atol=10**exponent)

    # fmt: on


class TestWildEdit:
    def test_wild_edit_pass(self, request):
        expected_dataset = idata.cnv_to_instrument_data(
            test_data / "19plus_V2_CTD-processing_example_wild_edit.cnv"
        )
        expected_conductivity = expected_dataset.measurements["c0S/m"].values

        dataset = idata.cnv_to_instrument_data(
            test_data / "19plus_V2_CTD-processing_example.cnv"
        )
        conductivity = dataset.measurements["c0S/m"].values
        flags = dataset.measurements["flag"].values

        wild_edit_output = p.wild_edit(conductivity, flags, 2, 20, 100, 0, False)

        request.node.return_value = wild_edit_output.tolist()
        assert np.all(wild_edit_output == expected_conductivity)


class TestWindowFilter:
    file_prefix = Path(test_data / "19plus_V2_CTD-processing_example")
    cnvdata = idata.cnv_to_instrument_data(f"{file_prefix}.cnv")
    pressure = cnvdata.measurements["prdM"].values
    flags = cnvdata.measurements["flag"].values
    window_width = 5
    half_width = 1  # only applies to gaussian
    offset = 0.25  # only applies to gaussian
    sample_interval = cnvdata.interval_s  # only applies to gaussian

    def test_boxcar_filter(self, request):
        expected_dataset = idata.cnv_to_instrument_data(f"{self.file_prefix}_boxcar_5.cnv")
        expected_pressure = expected_dataset.measurements["prdM"].values

        filtered_pressure = p.window_filter(
            self.pressure,
            self.flags,
            p.WindowFilterType.BOXCAR,
            self.window_width,
            self.cnvdata.interval_s,
            self.half_width,
            self.offset,
            False,
        )

        request.node.return_value = filtered_pressure.tolist()
        assert close_enough(filtered_pressure, expected_pressure, 3, 1e-12)

    def test_boxcar_filter_exclude_flags(self, request):
        expected_dataset = idata.cnv_to_instrument_data(f"{self.file_prefix}_boxcar_5_excluded.cnv")
        expected_pressure = expected_dataset.measurements["prdM"].values

        filtered_pressure = p.window_filter(
            self.pressure,
            self.flags,
            p.WindowFilterType.BOXCAR,
            self.window_width,
            self.cnvdata.interval_s,
            self.half_width,
            self.offset,
            True,
        )

        request.node.return_value = filtered_pressure.tolist()
        assert close_enough(filtered_pressure, expected_pressure, 3, 1e-12)

    def test_cosine_filter(self, request):
        expected_dataset = idata.cnv_to_instrument_data(f"{self.file_prefix}_cosine_5.cnv")
        expected_pressure = expected_dataset.measurements["prdM"].values

        filtered_pressure = p.window_filter(
            self.pressure,
            self.flags,
            p.WindowFilterType.COSINE,
            self.window_width,
            self.cnvdata.interval_s,
            self.half_width,
            self.offset,
            False,
        )

        request.node.return_value = filtered_pressure.tolist()
        assert close_enough(filtered_pressure, expected_pressure, 3, 1e-12)

    def test_cosine_filter_exclude_flags(self, request):
        expected_dataset = idata.cnv_to_instrument_data(f"{self.file_prefix}_cosine_5_excluded.cnv")
        expected_pressure = expected_dataset.measurements["prdM"].values

        filtered_pressure = p.window_filter(
            self.pressure,
            self.flags,
            p.WindowFilterType.COSINE,
            self.window_width,
            self.cnvdata.interval_s,
            self.half_width,
            self.offset,
            True,
        )

        request.node.return_value = filtered_pressure.tolist()
        assert close_enough(filtered_pressure, expected_pressure, 3, 1e-12)

    def test_triangle_filter(self, request):
        expected_dataset = idata.cnv_to_instrument_data(f"{self.file_prefix}_triangle_5.cnv")
        expected_pressure = expected_dataset.measurements["prdM"].values

        filtered_pressure = p.window_filter(
            self.pressure,
            self.flags,
            p.WindowFilterType.TRIANGLE,
            self.window_width,
            self.cnvdata.interval_s,
            self.half_width,
            self.offset,
            False,
        )

        request.node.return_value = filtered_pressure.tolist()
        assert close_enough(filtered_pressure, expected_pressure, 3, 1e-12)

    def test_triangle_filter_exclude_flags(self, request):
        expected_dataset = idata.cnv_to_instrument_data(f"{self.file_prefix}_triangle_5_excluded.cnv")
        expected_pressure = expected_dataset.measurements["prdM"].values

        filtered_pressure = p.window_filter(
            self.pressure,
            self.flags,
            p.WindowFilterType.TRIANGLE,
            self.window_width,
            self.cnvdata.interval_s,
            self.half_width,
            self.offset,
            True,
        )

        request.node.return_value = filtered_pressure.tolist()
        assert close_enough(filtered_pressure, expected_pressure, 3, 1e-12)

    def test_gaussian_filter(self, request):
        expected_dataset = idata.cnv_to_instrument_data(f"{self.file_prefix}_gaussian_5_1_025.cnv")
        expected_pressure = expected_dataset.measurements["prdM"].values

        filtered_pressure = p.window_filter(
            self.pressure,
            self.flags,
            p.WindowFilterType.GAUSSIAN,
            self.window_width,
            self.cnvdata.interval_s,
            self.half_width,
            self.offset,
            False,
        )

        request.node.return_value = filtered_pressure.tolist()
        assert close_enough(filtered_pressure, expected_pressure, 3, 1e-12)

    def test_gaussian_filter_exclude_flags(self, request):
        expected_dataset = idata.cnv_to_instrument_data(f"{self.file_prefix}_gaussian_5_1_025_excluded.cnv")
        expected_pressure = expected_dataset.measurements["prdM"].values

        filtered_pressure = p.window_filter(
            self.pressure,
            self.flags,
            p.WindowFilterType.GAUSSIAN,
            self.window_width,
            self.cnvdata.interval_s,
            self.half_width,
            self.offset,
            True,
        )

        request.node.return_value = filtered_pressure.tolist()
        assert close_enough(filtered_pressure, expected_pressure, 3, 1e-12)

    def test_median_filter(self, request):
        expected_dataset = idata.cnv_to_instrument_data(f"{self.file_prefix}_median_5.cnv")
        expected_pressure = expected_dataset.measurements["prdM"].values

        filtered_pressure = p.window_filter(
            self.pressure,
            self.flags,
            p.WindowFilterType.MEDIAN,
            self.window_width,
            self.cnvdata.interval_s,
            self.half_width,
            self.offset,
            False,
        )
        request.node.return_value = filtered_pressure.tolist()
        assert(close_enough(filtered_pressure, expected_pressure, 3, 1e-12))

    def test_median_filter_exclude_flags(self, request):
        expected_dataset = idata.cnv_to_instrument_data(f"{self.file_prefix}_median_5_excluded.cnv")
        expected_pressure = expected_dataset.measurements["prdM"].values

        filtered_pressure = p.window_filter(
            self.pressure,
            self.flags,
            p.WindowFilterType.MEDIAN,
            self.window_width,
            self.cnvdata.interval_s,
            self.half_width,
            self.offset,
            True,
        )
        request.node.return_value = filtered_pressure.tolist()
        assert(close_enough(filtered_pressure, expected_pressure, 3, 1e-12))


class TestBuoyancy:
    # fmt: off
    # Testing data comes from a CalCOFI cruise
    temperature = np.asarray(
        [16.7373, 16.5030, 16.1106, 14.3432, 13.0211, 12.0935, 11.3933, 11.2466, 10.9219, 10.4762, 9.9460]
    )
    pressure = np.asarray([10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110])
    salinity = np.asarray(
        [33.2410, 33.2321, 33.2091, 33.1329, 33.0762, 33.1391, 33.2560, 33.4015, 33.5683, 33.6766, 33.7794]
    )
    expected_N2_win30 = np.asarray(
        [-9.990e-29, 5.7702e-05, 1.9197e-04, 2.6735e-04, 2.1888e-04, 2.1620e-04, 1.7374e-04, 1.5761e-04, 1.6905e-04, 1.6099e-04, -9.990e-29]
    )
    expected_N_win30 = np.asarray(
        [-9.990e-29, 4.35, 7.94, 9.37, 8.48, 8.42, 7.55, 7.19, 7.45, 7.27, -9.990e-29]
    )
    expected_E_win30 = np.asarray(
        [-9.990e-29, 5.8901e-06, 1.9596e-05, 2.7290e-05, 2.2342e-05, 2.2069e-05, 1.7735e-05, 1.6089e-05, 1.7256e-05, 1.6433e-05, -9.990e-29]
    )
    expected_E_pow_8_win30 = np.asarray(
        [-9.990e-29, 589.0, 1959.6, 2729.0, 2234.2, 2206.9, 1773.5, 1608.9, 1725.6, 1643.3, -9.990e-29]
    )
    # fmt: on

    def test_buoyancy(self, request):
        output_dataframe = p.buoyancy(
            self.temperature,
            self.salinity,
            self.pressure,
            np.asarray([34.034167]),  # converted from metadata 34.02.03 N in H,M,S
            np.asarray([121.060556]),  # converted from metadata 121 03.38 W in H, M, S
            30,  # window size
            True,
        )

        request.node.return_value = {
            'N2': output_dataframe["N2"].to_list(),
            'N': output_dataframe["N"].to_list(),
            'E': output_dataframe["E"].to_list(),
            'E10^-8': output_dataframe["E10^-8"].to_list(),
        }
        # Comparing EOS-80 to TEOS-10 buoyancy calculations.
        # We do not expect them to agree better than +/-1.5% due to differences in the algorithms
        rel_tol = 0.015  # 1.5%
        assert output_dataframe["N2"].to_numpy() == pytest.approx(
            self.expected_N2_win30, rel=rel_tol
        )
        assert output_dataframe["N"].to_numpy() == pytest.approx(
            self.expected_N_win30, rel=rel_tol
        )
        assert output_dataframe["E"].to_numpy() == pytest.approx(
            self.expected_E_win30, rel=rel_tol
        )
        assert output_dataframe["E10^-8"].to_numpy() == pytest.approx(
            self.expected_E_pow_8_win30, rel=rel_tol
        )

    def test_buoyancy_eos80(self, request):
        output_dataframe = p.buoyancy(
            self.temperature,
            self.salinity,
            self.pressure,
            np.asarray([34.034167]),  # converted from metadata 34.02.03 N in H,M,S
            np.asarray([121.060556]),  # converted from metadata 121 03.38 W in H, M, S
            30,  # window size
            False,
        )

        request.node.return_value = {
            'N2': output_dataframe["N2"].to_list(),
            'N': output_dataframe["N"].to_list(),
            'E': output_dataframe["E"].to_list(),
            'E10^-8': output_dataframe["E10^-8"].to_list(),
        }
        # Comparing SBE Data Processing C++ to local Python results using the same EOS-80 calculations.
        # We expect very very close agreement: << 1% differnce
        rel_tol = 0.0026  # 0.26%
        assert output_dataframe["N2"].to_numpy() == pytest.approx(
            self.expected_N2_win30, rel=rel_tol
        )
        assert output_dataframe["N"].to_numpy() == pytest.approx(
            self.expected_N_win30, rel=rel_tol
        )
        assert output_dataframe["E"].to_numpy() == pytest.approx(
            self.expected_E_win30, rel=rel_tol
        )
        assert output_dataframe["E10^-8"].to_numpy() == pytest.approx(
            self.expected_E_pow_8_win30, rel=rel_tol
        )


class TestNitrate:
    dac_min = -5
    dac_max = 100
    voltages = np.array([0.095, 1, 2, 3, 4.095])

    def test_convert_nitrate_umno3(self, request):
        expected_umno3 = np.array([-5, 18.75625, 45.00625, 71.25625, 100])
        nitrate = c.convert_nitrate(self.voltages, dac_min=self.dac_min, dac_max=self.dac_max)
        request.node.return_value = nitrate.tolist()
        assert np.allclose(expected_umno3, nitrate, atol=0.000001)

    def test_convert_nitrate_mgnl(self, request):
        expected_mgnl = np.array([-0.070035, 0.26271879375, 0.63040254375, 0.99808629375, 1.4007])
        nitrate = c.convert_nitrate(self.voltages, dac_min=self.dac_min, dac_max=self.dac_max, units='mgNL')
        request.node.return_value = nitrate.tolist()
        assert np.allclose(expected_mgnl, nitrate, atol=0.000001)