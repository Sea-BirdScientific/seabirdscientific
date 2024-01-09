#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""TODO: countour docstring"""

# Native imports
from dataclasses import dataclass

# Third-party imports
import gsw
import numpy as np

# Sea-Bird imports

# Internal imports


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


def contour_from_t_s_p(temperature_C, salinity_PSU, pressure_dbar, min_salinity=0, lat=0, lon=0) -> ContourData:
    """Converts temperature T, salinity S, and pressure P to conservative temperature CT,
    absolute salinity SA, and potential density PD. CT is derived from ITS-90 temperature
    and practical salinity measurements. PD is derived from SA and CT.

    Args:
        temperature_C (ndarray): Measured temperature in degrees C
        salinity_PSU (ndarray): Measured salinity in practical salinity units
        pressure_dbar (ndarray): Measured pressure in decibar
        min_salinity (optional, float): Minimum salinity to include in contour data
        lat, lon (optional, float) : Used to determine absolute salinity (SA)

    Returns:
        ContourData: dataclass with xyz data for creating a TS plot
    """

    # Filter out data points where salinity < min_salinity
    salt_mask = salinity_PSU > min_salinity
    temperature_C = temperature_C[salt_mask]
    salinity_PSU = salinity_PSU[salt_mask]
    pressure_dbar = pressure_dbar[salt_mask]

    # Compute TEOS-10 quantities: SA, CT, potential_density
    absolute_salinity = gsw.SA_from_SP(salinity_PSU, pressure_dbar, lon, lat)
    conservative_temperature = gsw.CT_from_t(absolute_salinity, temperature_C, pressure_dbar)
    potential_density = gsw.rho(absolute_salinity, conservative_temperature, 0)  # TODO: parameterize reference density

    # Figure out T-S grid boundaries (mins and maxes)
    salt_min = absolute_salinity.min() - (0.01 * absolute_salinity.min())
    salt_max = absolute_salinity.max() + (0.01 * absolute_salinity.max())
    temperature_min = conservative_temperature.min() - (0.1 * conservative_temperature.min())
    temperature_max = conservative_temperature.max() + (0.1 * conservative_temperature.max())

    # Calculate how many grid cells we need in the x and y dimensions
    x_range = round((salt_max - salt_min) / 0.1 + 1, 0)
    y_range = round((temperature_max - temperature_min) + 1, 0)
    x_cells = x_range.astype(int)
    y_cells = y_range.astype(int)

    # Create conservative_temperature and absolute_salinity vectors of appropriate dimensions
    temperature_vector = np.linspace(1, y_range - 1, y_cells) + temperature_min
    salinity_vector = np.linspace(1, x_range - 1, x_cells) * 0.1 + salt_min

    # Loop to fill in density
    potential_density_matrix = np.zeros((y_cells, x_cells))
    for j in range(0, y_cells):
        potential_density_matrix[j, :] = gsw.rho(salinity_vector, temperature_vector[j], 0)

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


def contour_from_t_c_p(temperature_C, conductivity_mScm, pressure_dbar, min_salinity=0, lat=0, lon=0) -> ContourData:
    """Converts conductivity to salinity then calls derive_ct_sa_pd_from_t_s_p to
    derive conservative temperature, absolute salinity, and potential density.

    Args:
        temperature_C (ndarray): Measured temperature in degrees C
        Conductivity_mScm (ndarray): Measured conductivity in mSiemens/cm
        pressure_dbar (ndarray): Measured pressure in decibar
        min_salinity (optional, float): Minimum salinity to include in contour data
        lat, lon (optional, float) : Used to determine absolute salinity (SA)

    Returns:
        ContourData: dataclass with xyz data for creating a TS plot
    """
    salinity_PSU = gsw.SP_from_C(conductivity_mScm, temperature_C, pressure_dbar)
    return contour_from_t_s_p(temperature_C, salinity_PSU, pressure_dbar, min_salinity, lat, lon)
