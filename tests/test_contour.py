"""TS contour unit tests."""

import gsw
import numpy as np
import pytest

import seabirdscientific.contour as sc

# The CalCOFI cast used by TestBuoyancy in test_conversion.py
TEMPERATURE = np.array(
    [13.4288, 10.3983, 9.2891, 8.3246, 7.6800, 7.1504, 6.7090, 6.1575, 5.8453, 5.6158]
)
SALINITY = np.array(
    [33.1678, 33.6430, 33.9238, 33.9917, 34.0375, 34.0499, 34.0697, 34.1113, 34.1726, 34.2101]
)
PRESSURE = np.array([50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0])
LATITUDE = 34.034167
LONGITUDE = -121.060556


def expected_teos10(salinity, reference_pressure=0):
    """Absolute salinity, conservative temperature and potential density
    (as sigma) computed directly with gsw
    """
    absolute_salinity = gsw.SA_from_SP(salinity, PRESSURE, LONGITUDE, LATITUDE)
    conservative_temperature = gsw.CT_from_t(absolute_salinity, TEMPERATURE, PRESSURE)
    sigma = gsw.rho(absolute_salinity, conservative_temperature, reference_pressure) - 1000
    return absolute_salinity, conservative_temperature, sigma


class TestContourFromTSP:
    @pytest.mark.parametrize("reference_pressure", [0, 1000])
    def test_data_matches_teos10(self, reference_pressure):
        result = sc.contour_from_t_s_p(
            TEMPERATURE,
            SALINITY,
            PRESSURE,
            lat=LATITUDE,
            lon=LONGITUDE,
            reference_pressure=reference_pressure,
        )
        absolute_salinity, conservative_temperature, sigma = expected_teos10(
            SALINITY, reference_pressure
        )

        assert np.allclose(result.x, absolute_salinity, rtol=0, atol=1e-12)
        assert np.allclose(result.y, conservative_temperature, rtol=0, atol=1e-12)
        assert np.allclose(result.z, sigma, rtol=0, atol=1e-12)

    @pytest.mark.parametrize("reference_pressure", [0, 1000])
    def test_grid_is_density_of_grid_points(self, reference_pressure):
        result = sc.contour_from_t_s_p(
            TEMPERATURE,
            SALINITY,
            PRESSURE,
            lat=LATITUDE,
            lon=LONGITUDE,
            reference_pressure=reference_pressure,
        )

        assert result.z_mat.shape == (len(result.y_vec), len(result.x_vec))
        assert np.all(np.diff(result.x_vec) > 0)
        assert np.all(np.diff(result.y_vec) > 0)
        # rows are temperature, columns are salinity
        expected = (
            gsw.rho(result.x_vec[np.newaxis, :], result.y_vec[:, np.newaxis], reference_pressure)
            - 1000
        )
        assert np.allclose(result.z_mat, expected, rtol=0, atol=1e-12)

    # Coldest water above 10 C, below 10 C, and below 0 C. The old percentage-based padding
    # only covered the first
    @pytest.mark.parametrize(
        "temperature",
        [
            np.array([25.0, 20.0, 15.0, 12.0]),
            np.array([13.4, 9.3, 7.7, 5.6]),
            np.array([1.0, -0.5, -1.2, -1.8]),
        ],
        ids=["warm", "temperate", "polar"],
    )
    def test_grid_covers_data_with_fixed_padding(self, temperature):
        salinity = np.array([33.2, 33.9, 34.1, 34.3])
        pressure = np.array([10.0, 100.0, 300.0, 500.0])

        result = sc.contour_from_t_s_p(temperature, salinity, pressure)

        assert np.isclose(result.y_vec[0], np.nanmin(result.y) - sc.TEMPERATURE_GRID_PADDING)
        assert result.y_vec[-1] >= np.nanmax(result.y) + sc.TEMPERATURE_GRID_PADDING
        assert np.isclose(result.x_vec[0], np.nanmin(result.x) - sc.SALINITY_GRID_PADDING)
        assert result.x_vec[-1] >= np.nanmax(result.x) + sc.SALINITY_GRID_PADDING

    def test_grid_steps_are_exact(self):
        result = sc.contour_from_t_s_p(
            TEMPERATURE, SALINITY, PRESSURE, lat=LATITUDE, lon=LONGITUDE
        )

        assert np.allclose(np.diff(result.y_vec), sc.TEMPERATURE_GRID_STEP, rtol=0, atol=1e-12)
        assert np.allclose(np.diff(result.x_vec), sc.SALINITY_GRID_STEP, rtol=0, atol=1e-12)

    def test_min_salinity_excludes_samples(self):
        result = sc.contour_from_t_s_p(
            TEMPERATURE, SALINITY, PRESSURE, min_salinity=33.5, lat=LATITUDE, lon=LONGITUDE
        )
        excluded = SALINITY <= 33.5

        assert excluded.sum() == 1
        for values in (result.x, result.y, result.z):
            assert np.all(np.isnan(values[excluded]))
            assert np.all(np.isfinite(values[~excluded]))

        # excluded samples don't set the grid bounds
        full = sc.contour_from_t_s_p(TEMPERATURE, SALINITY, PRESSURE, lat=LATITUDE, lon=LONGITUDE)
        assert result.x_vec[0] > full.x_vec[0]

    def test_deprecated_keywords(self):
        expected = sc.contour_from_t_s_p(TEMPERATURE, SALINITY, PRESSURE)

        with pytest.warns(DeprecationWarning) as record:
            result = sc.contour_from_t_s_p(
                None,
                None,
                None,
                temperature_C=TEMPERATURE,
                salinity_PSU=SALINITY,
                pressure_dbar=PRESSURE,
            )

        messages = {str(warning.message) for warning in record}
        assert messages == {
            "Deprecated, use temperature",
            "Deprecated, use salinity",
            "Deprecated, use pressure",
        }

        assert np.array_equal(result.z_mat, expected.z_mat)
        assert np.array_equal(result.z, expected.z)


class TestContourFromTCP:
    # conductivity in mS/cm, from the practical salinity of the profile
    conductivity = gsw.C_from_SP(SALINITY, TEMPERATURE, PRESSURE)

    def test_matches_t_s_p_with_derived_salinity(self):
        result = sc.contour_from_t_c_p(
            TEMPERATURE,
            self.conductivity,
            PRESSURE,
            lat=LATITUDE,
            lon=LONGITUDE,
            reference_pressure=1000,
        )
        expected = sc.contour_from_t_s_p(
            TEMPERATURE,
            gsw.SP_from_C(self.conductivity, TEMPERATURE, PRESSURE),
            PRESSURE,
            lat=LATITUDE,
            lon=LONGITUDE,
            reference_pressure=1000,
        )

        for field in ("x", "y", "z", "x_vec", "y_vec", "z_mat"):
            assert np.array_equal(getattr(result, field), getattr(expected, field)), field

    def test_deprecated_keywords(self):
        expected = sc.contour_from_t_c_p(TEMPERATURE, self.conductivity, PRESSURE)

        with pytest.warns(DeprecationWarning) as record:
            result = sc.contour_from_t_c_p(
                None,
                None,
                None,
                temperature_C=TEMPERATURE,
                conductivity_mScm=self.conductivity,
                pressure_dbar=PRESSURE,
            )

        messages = {str(warning.message) for warning in record}
        assert messages == {
            "Deprecated, use temperature",
            "Deprecated, use conductivity",
            "Deprecated, use pressure",
        }

        assert np.array_equal(result.z_mat, expected.z_mat)


class TestGridVector:
    def test_exact_multiple_of_step_has_no_extra_point(self):
        # span of exactly 2.0 with padding 0 and step 0.5 is 5 points, 0.0 to 2.0
        result = sc._grid_vector(np.array([0.0, 2.0]), step=0.5, padding=0)

        assert np.allclose(result, [0.0, 0.5, 1.0, 1.5, 2.0], rtol=0, atol=1e-12)

    def test_reaches_past_maximum_and_ignores_nan(self):
        result = sc._grid_vector(np.array([np.nan, 0.0, 1.05]), step=0.5, padding=0.1)

        assert np.isclose(result[0], -0.1)
        assert result[-1] >= 1.15
        assert result[-2] < 1.15
