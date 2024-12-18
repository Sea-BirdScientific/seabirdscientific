#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""A collection of calibration data from SBS instruments.

These are calibration coefficients from specific SBS intstruments that are primarily used for test purposes.
"""

# Native imports
from dataclasses import dataclass

# Third-party imports

# Sea-Bird imports

# Internal imports


class TemperatureCoefficients:
    """
    Args:
        a0 (float): a0 calibration coefficient for the temperature sensor
        a1 (float): a1 calibration coefficient for the temperature sensor
        a2 (float): a2 calibration coefficient for the temperature sensor
        a3 (float): a3 calibration coefficient for the temperature sensor
    """

    def __init__(self, a0, a1, a2, a3):
        self.a0 = a0
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3


class PressureCoefficients:
    """
    Args:
        pa0 (float): pa0 calibration coefficient for the pressure sensor
        pa1 (float): pa1 calibration coefficient for the pressure sensor
        pa2 (float): pa2 calibration coefficient for the pressure sensor
        ptca0 (float): ptca0 calibration coefficient for the pressure sensor
        ptca1 (float): ptca1 calibration coefficient for the pressure sensor
        ptca2 (float): ptca2 calibration coefficient for the pressure sensor
        ptcb0 (float): ptcb0 calibration coefficient for the pressure sensor
        ptcb1 (float): ptcb1 calibration coefficient for the pressure sensor
        ptcb2 (float): ptcb2 calibration coefficient for the pressure sensor
        ptempa0 (float): ptempa0 calibration coefficient for the pressure sensor
        ptempa1 (float): ptempa1 calibration coefficient for the pressure sensor
        ptempa2 (float): ptempa2 calibration coefficient for the pressure sensor
    """

    def __init__(self, pa0, pa1, pa2, ptca0, ptca1, ptca2, ptcb0, ptcb1, ptcb2, ptempa0, ptempa1, ptempa2):
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
    Args:
        g (float): g calibration coefficient for the conductivity sensor
        h (float): h calibration coefficient for the conductivity sensor
        i (float): i calibration coefficient for the conductivity sensor
        j (float): j calibration coefficient for the conductivity sensor
        cpcor (float): cpcor calibration coefficient for the conductivity sensor
        ctcor (float): ctcor calibration coefficient for the conductivity sensor
        wbotc (float): Wien bridge oscillator temperature coefficient
            see the 37 Manual: https://www.seabird.com/asset-get.download.jsa?id=54627862348
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
    Args:
        slope (float): units/count for digital, units/V for analog
        offset (float): dark counts for digital, dark voltage for analog
    """

    def __init__(self, slope, offset):
        self.slope = slope
        self.offset = offset


class Oxygen63Coefficients:
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
    Args:
        ta0 (float): calibration coefficient for the thermistor in the SBE63 sensor
        ta1 (float): calibration coefficient for the thermistor in the SBE63 sensor
        ta2 (float): calibration coefficient for the thermistor in the SBE63 sensor
        ta3 (float): calibration coefficient for the thermistor in the SBE63 sensor
    """

    def __init__(self, ta0, ta1, ta2, ta3):
        self.ta0 = ta0
        self.ta1 = ta1
        self.ta2 = ta2
        self.ta3 = ta3


class Oxygen43Coefficients:
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
    Args:
        offset (float): calibration offset
        slope (float): calibration slope
    """

    def __init__(self, slope, offset):
        self.slope = slope
        self.offset = offset


class PARCoefficients:
    """
    Args:
        im (float): immersion coefficient
        a0 (float): calibration slope
        a1 (float): calibration offset
        multiplier (float): 1.0 for units of μEinsteins/m2 sec
    """

    def __init__(self, im, a0, a1, multiplier):
        self.im = im
        self.a0 = a0
        self.a1 = a1
        self.multiplier = multiplier
