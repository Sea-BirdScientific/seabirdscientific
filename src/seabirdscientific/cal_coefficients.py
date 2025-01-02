#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""Calibration coefficient objects for data conversion"""

# Native imports

# Third-party imports

# Sea-Bird imports

# Internal imports


class TemperatureCoefficients:
    """
    :param a0: calibration coefficient for the temperature sensor
    :param a1: calibration coefficient for the temperature sensor
    :param a2: calibration coefficient for the temperature sensor
    :param a3: calibration coefficient for the temperature sensor
    """

    def __init__(self, a0, a1, a2, a3):
        self.a0 = a0
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3


class PressureCoefficients:
    """
    :param pa0: calibration coefficient for the pressure sensor
    :param pa1: calibration coefficient for the pressure sensor
    :param pa2: calibration coefficient for the pressure sensor
    :param ptca0: calibration coefficient for the pressure sensor
    :param ptca1: calibration coefficient for the pressure sensor
    :param ptca2: calibration coefficient for the pressure sensor
    :param ptcb0: calibration coefficient for the pressure sensor
    :param ptcb1: calibration coefficient for the pressure sensor
    :param ptcb2: calibration coefficient for the pressure sensor
    :param ptempa0: calibration coefficient for the pressure sensor
    :param ptempa1: calibration coefficient for the pressure sensor
    :param ptempa2: calibration coefficient for the pressure sensor
    """

    def __init__(
        self, pa0, pa1, pa2, ptca0, ptca1, ptca2, ptcb0, ptcb1, ptcb2, ptempa0, ptempa1, ptempa2
    ):
        self.pa0 = pa0
        self.pa1 = pa1
        self.pa2 = pa2
        self.ptca0 = ptca0
        self.ptca1 = ptca1
        self.ptca2 = ptca2
        self.ptcb0 = ptcb0
        self.ptcb1 = ptcb1
        self.ptcb2 = ptcb2
        self.ptempa0 = ptempa0
        self.ptempa1 = ptempa1
        self.ptempa2 = ptempa2


class ConductivityCoefficients:
    """
    :param g: calibration coefficient for the conductivity sensor
    :param h: calibration coefficient for the conductivity sensor
    :param i: calibration coefficient for the conductivity sensor
    :param j: calibration coefficient for the conductivity sensor
    :param cpcor: calibration coefficient for the conductivity sensor
    :param ctcor: calibration coefficient for the conductivity sensor
    :param wbotc: bridge oscillator temperature coefficient see the
        37 Manual: https://www.seabird.com/asset-get.download.jsa?id=54627862348
    """
    def __init__(self, g, h, i, j, cpcor, ctcor, wbotc):
        self.g = g
        self.h = h
        self.i = i
        self.j = j
        self.cpcor = cpcor
        self.ctcor = ctcor
        self.wbotc = wbotc


class ECOCoefficients:
    """
    :param slope: units/count for digital, units/V for analog
    :param offset: dark counts for digital, dark voltage for analog
    """
    def __init__(self, slope, offset):
        self.slope = slope
        self.offset = offset


class Oxygen63Coefficients:
    """
    :param a0: calibration coefficient
    :param a1: calibration coefficient
    :param a2: calibration coefficient
    :param b0: calibration coefficient
    :param b1: calibration coefficient
    :param c0: calibration coefficient
    :param c1: calibration coefficient
    :param c2: calibration coefficient
    :param e: calibration coefficient
    """
    def __init__(self, a0, a1, a2, b0, b1, c0, c1, c2, e):
        self.a0 = a0
        self.a1 = a1
        self.a2 = a2
        self.b0 = b0
        self.b1 = b1
        self.c0 = c0
        self.c1 = c1
        self.c2 = c2
        self.e = e


class Thermistor63Coefficients:
    """
    :param ta0: calibration coefficient for the thermistor in the SBE63 sensor
    :param ta1: calibration coefficient for the thermistor in the SBE63 sensor
    :param ta2: calibration coefficient for the thermistor in the SBE63 sensor
    :param ta3: calibration coefficient for the thermistor in the SBE63 sensor
    """
    def __init__(self, ta0, ta1, ta2, ta3):
        self.ta0 = ta0
        self.ta1 = ta1
        self.ta2 = ta2
        self.ta3 = ta3


class Oxygen43Coefficients:
    """
    :param soc: linear scaling calibration coefficient
    :param v_offset: voltage at zero oxygen signal
    :param tau_20: sensor time constant tau(T,P) at 20 C, 1 atmosphere, 0 PSU;
        slope term in calculation of tau(T,P)
    :param a: calibration coefficient
    :param b: calibration coefficient
    :param c: calibration coefficient
    :param e: calibration coefficient
    :param d0: calibration terms used in calculation of tau(T,P)
    :param d1: calibration terms used in calculation of tau(T,P)
    :param d2: calibration terms used in calculation of tau(T,P)
    :param h1: calibration terms used for hysteresis correction
    :param h2: calibration terms used for hysteresis correction
    :param h3: calibration terms used for hysteresis correction
    """
    def __init__(self, soc, v_offset, tau_20, a, b, c, e, d0, d1, d2, h1, h2, h3):
        self.soc = soc
        self.v_offset = v_offset
        self.tau_20 = tau_20
        self.a = a
        self.b = b
        self.c = c
        self.e = e
        self.d0 = d0  # TODO: the sensor lists a d0 coefficient, but it doesn't seem to be used?
        self.d1 = d1
        self.d2 = d2
        self.h1 = h1
        self.h2 = h2
        self.h3 = h3


class PH18Coefficients:
    """
    :param offset: calibration offset
    :param slope: calibration slope
    """

    def __init__(self, slope, offset):
        self.slope = slope
        self.offset = offset


class PARCoefficients:
    """
    :param im: immersion coefficient
    :param a0: calibration slope
    :param a1: calibration offset
    :param multiplier: 1.0 for units of μEinsteins/m2 sec
    """
    def __init__(self, im, a0, a1, multiplier):
        self.im = im
        self.a0 = a0
        self.a1 = a1
        self.multiplier = multiplier


class PHSeaFETInternalCoefficients:
    """
    :param int_k0: Internal K0 calibration coefficient
    :param int_k2: Internal K2 calibration coefficient
    """

    def __init__(self, int_k0, int_k2):
        self.int_k0 = int_k0
        self.int_k2 = int_k2


class PHSeaFETExternalCoefficients:
    """
    :param ext_k0: External K0 calibration coefficient
    :param ext_k2: External K2 calibration coefficient
    :param f1: f(P) calibration coefficient
    :param f2: f(P) calibration coefficient
    :param f3: f(P) calibration coefficient
    :param f4: f(P) calibration coefficient
    :param f5: f(P) calibration coefficient
    :param f6: f(P) calibration coefficient
    """

    def __init__(self, ext_k0, ext_k2, f1, f2, f3, f4, f5, f6):
        self.ext_k0 = ext_k0
        self.ext_k2 = ext_k2
        self.f1 = f1
        self.f2 = f2
        self.f3 = f3
        self.f4 = f4
        self.f5 = f5
        self.f6 = f6
