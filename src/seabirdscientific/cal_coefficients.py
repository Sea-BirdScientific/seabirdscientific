#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""Calibration coefficient objects for data conversion"""

# Native imports
from dataclasses import dataclass
import warnings

# Third-party imports

# Sea-Bird imports

# Internal imports

@dataclass
class TemperatureCoefficients:
    """
    :param a0: calibration coefficient
    :param a1: calibration coefficient
    :param a2: calibration coefficient
    :param a3: calibration coefficient
    """
    a0: float
    a1: float
    a2: float
    a3: float


@dataclass
class PressureCoefficients:
    """
    :param pa0: calibration coefficient
    :param pa1: calibration coefficient
    :param pa2: calibration coefficient
    :param ptca0: calibration coefficient
    :param ptca1: calibration coefficient
    :param ptca2: calibration coefficient
    :param ptcb0: calibration coefficient
    :param ptcb1: calibration coefficient
    :param ptcb2: calibration coefficient
    :param ptempa0: calibration coefficient
    :param ptempa1: calibration coefficient
    :param ptempa2: calibration coefficient
    """
    pa0: float
    pa1: float
    pa2: float
    ptca0: float
    ptca1: float
    ptca2: float
    ptcb0: float
    ptcb1: float
    ptcb2: float
    ptempa0: float
    ptempa1: float
    ptempa2: float


@dataclass
class ConductivityCoefficients:
    """
    :param g: calibration coefficient
    :param h: calibration coefficient
    :param i: calibration coefficient
    :param j: calibration coefficient
    :param cpcor: compressibility coefficient
    :param ctcor: calibration coefficient
    :param wbotc: bridge oscillator temperature coefficient see the
        37 Manual: https://www.seabird.com/asset-get.download.jsa?id=54627862348
    """
    g: float
    h: float
    i: float
    j: float
    cpcor: float
    ctcor: float
    wbotc: float


@dataclass
class ECOCoefficients:
    """
    :param slope: units/count for digital, units/V for analog
    :param offset: dark counts for digital, dark voltage for analog
    """
    slope: float
    offset: float


@dataclass
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
    a0: float
    a1: float
    a2: float
    b0: float
    b1: float
    c0: float
    c1: float
    c2: float
    e: float


@dataclass
class Thermistor63Coefficients:
    """
    :param ta0: calibration coefficient
    :param ta1: calibration coefficient
    :param ta2: calibration coefficient
    :param ta3: calibration coefficient
    """
    ta0: float
    ta1: float
    ta2: float
    ta3: float


@dataclass
class Oxygen43Coefficients:
    """
    :param soc: linear scaling coefficient
    :param v_offset: voltage at zero oxygen signal
    :param tau_20: sensor time constant tau(T,P) at 20 C, 1 atmosphere, 0 PSU;
        slope term in calculation of tau(T,P)
    :param a: calibration coefficient
    :param b: calibration coefficient
    :param c: calibration coefficient
    :param e: calibration coefficient
    :param d0: calibration terms used to calculate tau(T,P)
    :param d1: calibration terms used to calculate tau(T,P)
    :param d2: calibration terms used to calculate tau(T,P)
    :param h1: calibration terms used for hysteresis correction
    :param h2: calibration terms used for hysteresis correction
    :param h3: calibration terms used for hysteresis correction
    """
    soc: float
    v_offset: float
    tau_20: float
    a: float
    b: float
    c: float
    e: float
    d0: float  # TODO: the sensor lists a d0 coefficient, but it doesn't seem to be used?
    d1: float
    d2: float
    h1: float
    h2: float
    h3: float


@dataclass
class PH18Coefficients:
    """
    :param slope: calibration slope
    :param offset: calibration offset
    """
    slope: float
    offset: float


@dataclass
class PARCoefficients:
    """
    :param im: immersion coefficient
    :param a0: calibration slope
    :param a1: calibration offset
    :param multiplier: 1.0 for units of μEinsteins/m2 sec
    """
    im: float
    a0: float
    a1: float
    multiplier: float


class PHSeaFETInternalCoefficients:
    """
    :param k0: K0 coefficient
    :param k2: K2 coefficient
    :param int_k0: Deprecated, use k0
    :param int_k2: Deprecated, use k2
    """
    def __init__(
            self,
            k0: float = None,
            k2: float = None,
            int_k0: float = None,
            int_k2: float = None,
    ):
        self.k0 = k0
        self.k2 = k2
        if int_k0 is not None:
            self.int_k0 = int_k0
            warnings.warn("int_k0 is deprecated, use k0", DeprecationWarning)
        if int_k2 is not None:
            self.int_k2 = int_k2
            warnings.warn("int_k2 is deprecated, use k2", DeprecationWarning)


class PHSeaFETExternalCoefficients:
    """
    :param k0: K0 coefficient
    :param k2: K2 coefficient
    :param ext_k0: Deprecated, use k0
    :param ext_k2: Deprecated, use k2
    :param f1: f(P) coefficient
    :param f2: f(P) coefficient
    :param f3: f(P) coefficient
    :param f4: f(P) coefficient
    :param f5: f(P) coefficient
    :param f6: f(P) coefficient
    """
    def __init__(
            self,
            k0: float = None,
            k2: float = None,
            ext_k0: float = None,
            ext_k2: float = None,
            f1: float = None,
            f2: float = None,
            f3: float = None,
            f4: float = None,
            f5: float = None,
            f6: float = None,
    ):
        self.k0 = k0
        self.k2 = k2
        self.f1 = f1
        self.f2 = f2
        self.f3 = f3
        self.f4 = f4
        self.f5 = f5
        self.f6 = f6

        if ext_k0 is not None:
            self.ext_k0 = ext_k0
            warnings.warn("int_k0 is deprecated, use k0", DeprecationWarning)
        if ext_k2 is not None:
            self.ext_k2 = ext_k2
            warnings.warn("int_k2 is deprecated, use k2", DeprecationWarning)
