"""Contains a class and functions to support calculation of temperature
and salinity (TS) contours.
"""

import warnings
from dataclasses import dataclass

import gsw
import numpy as np


@dataclass
class ContourData:
    """Container for contour data for a TS plot.

    Example variables are shown next to each field
    """

    x: np.ndarray  # absolute salinity
    y: np.ndarray  # conservative temperature
    z: np.ndarray  # potential density
    x_vec: np.ndarray  # absolute salinity vector
    y_vec: np.ndarray  # conservative temperature vector
    z_mat: np.ndarray  # potential density matrix


def contour_from_t_s_p(
    temperature: np.ndarray,
    salinity: np.ndarray,
    pressure: np.ndarray,
    min_salinity: float = 0,
    lat: float = 0,
    lon: float = 0,
    reference_pressure: float = 0,
    temperature_C: np.ndarray = None,
    salinity_PSU: np.ndarray = None,
    pressure_dbar: np.ndarray = None,
) -> ContourData:
    """Converts temperature (T), salinity (S), and pressure (P) to
    conservative temperature (CT), absolute salinity (SA), and potential
    density (PD). CT is derived from ITS-90 temperature and practical
    salinity measurements. PD is derived from SA and CT.

    :param temperature_C: Measured temperature in degrees C
    :param salinity_PSU: Measured salinity in practical salinity units
    :param pressure_dbar: Measured pressure in decibars
    :param min_salinity: Minimum salinity to include in contour data.
        Defaults to 0
    :param lat: Used to determine absolute salinity (SA). Defaults to 0
    :param lon: Used to determine absolute salinity (SA). Defaults to 0

    :return: dataclass with xyz data for creating a TS plot
    """

    if temperature_C is not None:
        warnings.warn("Deprecated, use temperature", DeprecationWarning)
        temperature = temperature_C

    if salinity_PSU is not None:
        warnings.warn("Deprecated, use salinity", DeprecationWarning)
        salinity = salinity_PSU

    if pressure_dbar is not None:
        warnings.warn("Deprecated, use pressure", DeprecationWarning)
        pressure = pressure_dbar

    # Mark data as nan where salinity < min_salinity
    salinity_mask = salinity > min_salinity
    temperature = np.where(salinity_mask, temperature, np.nan)
    salinity = np.where(salinity_mask, salinity, np.nan)
    pressure = np.where(salinity_mask, pressure, np.nan)

    # Compute TEOS-10 quantities: SA, CT, potential_density
    absolute_salinity = gsw.SA_from_SP(salinity, pressure, lon, lat)
    conservative_temperature = gsw.CT_from_t(absolute_salinity, temperature, pressure)
    potential_density = gsw.rho(absolute_salinity, conservative_temperature, reference_pressure)

    # Figure out T-S grid boundaries (mins and maxes)
    min_s = np.nanmin(absolute_salinity) - (0.01 * np.nanmin(absolute_salinity))
    max_s = np.nanmax(absolute_salinity) + (0.01 * np.nanmax(absolute_salinity))
    min_t = np.nanmin(conservative_temperature) - (0.1 * np.nanmin(conservative_temperature))
    max_t = np.nanmax(conservative_temperature) + (0.1 * np.nanmax(conservative_temperature))

    # Calculate how many grid cells we need in the x and y dimensions
    x_range = round((max_s - min_s) / 0.1 + 1, 0)
    y_range = round((max_t - min_t) + 1, 0)
    x_cells = x_range.astype(int)
    y_cells = y_range.astype(int)

    # Create conservative_temperature and absolute_salinity vectors of appropriate dimensions
    temperature_vector = np.linspace(1, y_range - 1, y_cells) + min_t
    salinity_vector = np.linspace(1, x_range - 1, x_cells) * 0.1 + min_s

    # Loop to fill in density
    potential_density_matrix = np.zeros((y_cells, x_cells))
    for j in range(y_cells):
        potential_density_matrix[j, :] = gsw.rho(
            salinity_vector, temperature_vector[j], reference_pressure
        )

    # Subtract 1000 to convert to sigma-t
    potential_density_matrix -= 1000
    potential_density -= 1000

    contour_data = ContourData(
        x=absolute_salinity,
        y=conservative_temperature,
        z=potential_density,
        x_vec=salinity_vector,
        y_vec=temperature_vector,
        z_mat=potential_density_matrix,
    )

    return contour_data


def contour_from_t_c_p(
    temperature: np.ndarray,
    conductivity: np.ndarray,
    pressure: np.ndarray,
    min_salinity: float = 0,
    lat: float = 0,
    lon: float = 0,
    reference_pressure: float = 0,
    temperature_C: np.ndarray = None,
    conductivity_mScm: np.ndarray = None,
    pressure_dbar: np.ndarray = None,
) -> ContourData:
    """Converts conductivity (C) to salinity (S) then calls
    derive_ct_sa_pd_from_t_s_p to derive conservative temperature (CT),
    absolute salinity (SA), and potential density (PD).

    :param temperature_C: Measured temperature in degrees C
    :param Conductivity_mScm: Measured conductivity in mSiemens/cm
    :param pressure_dbar: Measured pressure in decibar
    :param min_salinity: Minimum salinity to include in contour data.
        Defaults to 0
    :param lat: Used to determine absolute salinity (SA). Defaults to 0
    :param lon: Used to determine absolute salinity (SA). Defaults to 0

    :return: dataclass with xyz data for creating a TS plot
    """

    if temperature_C is not None:
        warnings.warn("Deprecated, use temperature", DeprecationWarning)
        temperature = temperature_C

    if conductivity_mScm is not None:
        warnings.warn("Deprecated, use conductivity", DeprecationWarning)
        conductivity = conductivity_mScm

    if pressure_dbar is not None:
        warnings.warn("Deprecated, use pressure", DeprecationWarning)
        pressure = pressure_dbar

    salinity = gsw.SP_from_C(conductivity, temperature, pressure)
    return contour_from_t_s_p(
        temperature, salinity, pressure, min_salinity, lat, lon, reference_pressure
    )
