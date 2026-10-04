"""Deprecated EOS-80 processing function unit tests."""

import numpy as np
import pytest
import seawater as sw

import seabirdscientific.eos80_conversion as ec
import seabirdscientific.eos80_processing as ep

# A short profile from the CalCOFI cast used by TestBuoyancy in test_conversion.py
SALINITY = np.array([33.17, 33.64, 33.92, 34.04, 34.21])
TEMPERATURE = np.array([13.43, 10.40, 9.29, 7.68, 5.62])
PRESSURE = np.array([50.0, 100.0, 150.0, 250.0, 500.0])


class TestBouyancyFrequency:
    def test_delegates_to_eos80_conversion(self):
        gravity = 9.7963
        expected = ec.bouyancy_frequency(TEMPERATURE, SALINITY, PRESSURE, gravity)

        with pytest.warns(DeprecationWarning, match="eos80_conversion.buoyancy_frequency"):
            result = ep.bouyancy_frequency(TEMPERATURE, SALINITY, PRESSURE, gravity)

        assert result == expected


class TestDensity:
    def test_returns_seawater_density_as_sigma(self):
        # SBE Data Processing reports density as sigma (density - 1000 kg/m^3)
        expected = sw.dens(SALINITY, TEMPERATURE, PRESSURE) - 1000

        with pytest.warns(DeprecationWarning, match="seawater.dens"):
            result = ep.density(SALINITY, TEMPERATURE, PRESSURE)

        assert np.allclose(result, expected, rtol=0, atol=1e-12)

    def test_scalar_input_returns_array(self):
        with pytest.warns(DeprecationWarning):
            result = ep.density(35.0, 10.0, 0.0)

        assert result.shape == (1,)


class TestPotentialTemperature:
    # Expected values from the v2.8.1 implementation (a port of SeaSoft's PoTemp), which
    # integrates from in-situ pressure p0 to reference pressure pr
    @pytest.mark.parametrize(
        "reference_pressure, expected",
        [
            (0.0, [13.423129, 10.388356, 9.273658, 7.655573, 5.578004]),
            (1000.0, [13.566245, 10.510791, 9.388481, 7.75875, 5.66589]),
        ],
    )
    def test_matches_v2(self, reference_pressure, expected):
        with pytest.warns(DeprecationWarning, match="seawater.ptemp"):
            result = ep.potential_temperature(
                SALINITY, TEMPERATURE, PRESSURE, np.full(5, reference_pressure)
            )

        # seawater differs from the v2 port by ~1e-5 deg C
        assert np.allclose(result, expected, rtol=0, atol=1e-4)

    def test_uses_in_situ_pressure(self):
        # At its own reference pressure, potential temperature equals in-situ temperature
        with pytest.warns(DeprecationWarning):
            result = ep.potential_temperature(SALINITY, TEMPERATURE, PRESSURE, PRESSURE)

        assert np.allclose(result, TEMPERATURE, rtol=0, atol=1e-10)


class TestAdiabaticTemperatureGradient:
    def test_delegates_to_eos80_conversion(self):
        expected = ec.adiabatic_temperature_gradient(SALINITY, TEMPERATURE, PRESSURE)

        with pytest.warns(
            DeprecationWarning, match="eos80_conversion.adiabatic_temperature_gradient"
        ):
            result = ep.adiabatic_temperature_gradient(SALINITY, TEMPERATURE, PRESSURE)

        assert np.array_equal(result, expected)
