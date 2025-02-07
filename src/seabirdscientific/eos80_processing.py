#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""EOS80 functions to support legacy processing"""


import numpy as np
from scipy import stats


def bouyancy_frequency(
    temp_ITS_subset: np.ndarray,  # pylint: disable=invalid-name #TODO: change this to be all lower case for TKIT-75
    salinity_prac_subset: np.ndarray,
    pressure_dbar_subset: np.ndarray,
    gravity: float,
):
    """Calculates an N^2 value (buoyancy frequency) for the given window
    of temperature, salinity, and pressure, at the given latitude.

    Expects temperature as ITS-90 temperature, salinity as practical
    salinity, and pressure as dbar, all of the same length. Performs the
    calculation following the SBE Data Processing formula using E0S-80
    calculations for potential temp and density

    :param temp_ITS_subset: ITS-90 temperature values for the given
        window
    :param salinity_prac_subset: practical salinity values for the given
        window
    :param pressure_dbar_subset: pressure values for the given window
    :param gravity: gravity value

    :return: A single N^2 [Brunt-Väisälä (buoyancy) frequency]
    """

    db_to_pa = 1e4

    # Wrap these as a length-1 array so that GSW accepts them
    pressure_bar = np.array([np.mean(pressure_dbar_subset)])
    temperature_bar = np.array([np.mean(temp_ITS_subset)])
    salinity_bar = np.array([np.mean(salinity_prac_subset)])

    # Compute average density over the window
    # rho_bar0 = gsw.rho(salinity_bar, temperature_bar, pressure_bar)[0]
    rho_bar = density(salinity_bar, temperature_bar, pressure_bar)[0]

    # Use SBE DP (EOS-80) formulas for potential temp and density
    theta = potential_temperature(
        salinity_prac_subset, temp_ITS_subset, pressure_dbar_subset, pressure_bar
    )
    v_vals = 1.0 / density(salinity_prac_subset, theta, pressure_bar)

    # Estimate vertical gradient of specific volume
    dvdp_result = stats.linregress(pressure_dbar_subset, v_vals)

    # Compute EOS-80 N2 combining computed average density and vertical gradient
    # we index into v_bar, alpha_bar, and beta_bar as they are all arrays of len 1
    n2 = 0 - (rho_bar**2 * gravity**2 * dvdp_result.slope / db_to_pa)
    return n2


# pylint: disable=too-many-statements #TODO: break this function up
def density(s0: np.ndarray, t: np.ndarray, p0: np.ndarray) -> np.ndarray:
    """EOS-80 density calculation.

    This was ported from CSharedCalc::Density()

    :param s0: salinity data
    :param t: temperature data
    :param p0: pressure data

    :return: resulting density data
    """

    b0 = 8.24493e-1
    b1 = -4.0899e-3
    b2 = 7.6438e-5
    b3 = -8.2467e-7
    b4 = 5.3875e-9

    c0 = -5.72466e-3
    c1 = 1.0227e-4
    c2 = -1.6546e-6

    d0 = 4.8314e-4

    a0 = 999.842594
    a1 = 6.793952e-2
    a2 = -9.095290e-3
    a3 = 1.001685e-4
    a4 = -1.120083e-6
    a5 = 6.536332e-9

    fq0 = 54.6746
    fq1 = -0.603459
    fq2 = 1.09987e-2
    fq3 = -6.1670e-5

    g0 = 7.944e-2
    g1 = 1.6483e-2
    g2 = -5.3009e-4

    i0 = 2.2838e-3
    i1 = -1.0981e-5
    i2 = -1.6078e-6

    j0 = 1.91075e-4

    m0 = -9.9348e-7
    m1 = 2.0816e-8
    m2 = 9.1697e-10

    e0 = 19652.21
    e1 = 148.4206
    e2 = -2.327105
    e3 = 1.360477e-2
    e4 = -5.155288e-5

    h0 = 3.239908
    h1 = 1.43713e-3
    h2 = 1.16092e-4
    h3 = -5.77905e-7

    k0 = 8.50935e-5
    k1 = -6.12293e-6
    k2 = 5.2787e-8

    s0, t, p0 = np.broadcast_arrays(s0, t, p0)
    p = p0.copy()
    s = s0.copy()

    t2 = t * t
    t3 = t * t2
    t4 = t * t3
    t5 = t * t4
    s[s <= 0.0] = 0.000001
    s32 = s**1.5
    p /= 10.0
    sigma = (
        a0
        + a1 * t
        + a2 * t2
        + a3 * t3
        + a4 * t4
        + a5 * t5
        + (b0 + b1 * t + b2 * t2 + b3 * t3 + b4 * t4) * s
        + (c0 + c1 * t + c2 * t2) * s32
        + d0 * s * s
    )

    kw = e0 + e1 * t + e2 * t2 + e3 * t3 + e4 * t4
    aw = h0 + h1 * t + h2 * t2 + h3 * t3
    bw = k0 + k1 * t + k2 * t2

    k = (
        kw
        + (fq0 + fq1 * t + fq2 * t2 + fq3 * t3) * s
        + (g0 + g1 * t + g2 * t2) * s32
        + (aw + (i0 + i1 * t + i2 * t2) * s + (j0 * s32)) * p
        + (bw + (m0 + m1 * t + m2 * t2) * s) * p * p
    )

    val = 1 - p / k

    val[val <= 0] = np.nan
    val = sigma / val

    return val


def potential_temperature(
    s: np.ndarray, t0: np.ndarray, p0: np.ndarray, pr: np.ndarray
) -> np.ndarray:
    """EOS-80 potential temperature calculation.

    This was ported from CSharedCalc::PoTemp()

    :param s: sainity data
    :param t0: temperature data
    :param p0: subset pressure data
    :param pr: pressure data

    :return: calculated potential temperature data
    """

    s, t0, p0, pr = np.broadcast_arrays(s, t0, p0, pr)

    p = p0.copy()
    t = t0.copy()
    h = pr - p
    xk = h * adiabatic_temperature_gradient(s, t, p)
    t += 0.5 * xk
    q = xk
    p += 0.5 * h
    xk = h * adiabatic_temperature_gradient(s, t, p)
    t += 0.29289322 * (xk - q)
    q = 0.58578644 * xk + 0.121320344 * q
    xk = h * adiabatic_temperature_gradient(s, t, p)
    t += 1.707106781 * (xk - q)
    q = 3.414213562 * xk - 4.121320344 * q
    p += 0.5 * h
    xk = h * adiabatic_temperature_gradient(s, t, p)
    temp = t + (xk - 2.0 * q) / 6.0

    return temp


def adiabatic_temperature_gradient(
    s: np.ndarray,
    t: np.ndarray,
    p: np.ndarray,
) -> np.ndarray:
    """EOS-80 adiabatic lapse rate calculation.

    This was ported from CSharedCalc::ATG()

    :param s: salinity data
    :param t: temperature data
    :param p: pressure data

    :return: the resulting adiabatic lapse rate
    """

    s, t, p = np.broadcast_arrays(s, t, p)

    ds = s - 35.0
    atg = (
        (
            ((-2.1687e-16 * t + 1.8676e-14) * t - 4.6206e-13) * p
            + (
                (2.7759e-12 * t - 1.1351e-10) * ds
                + ((-5.4481e-14 * t + 8.733e-12) * t - 6.7795e-10) * t
                + 1.8741e-8
            )
        )
        * p
        + (-4.2393e-8 * t + 1.8932e-6) * ds
        + ((6.6228e-10 * t - 6.836e-8) * t + 8.5258e-6) * t
        + 3.5803e-5
    )

    return atg
