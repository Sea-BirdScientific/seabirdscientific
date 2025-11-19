"""Calibration coefficient dataclasses for data conversion."""

# Native imports
from dataclasses import dataclass
from typing import Optional
import warnings

# Third-party imports

# Sea-Bird imports

# Internal imports


@dataclass
class TemperatureCoefficients:
    """
    :param a0: coefficient
    :param a1: coefficient
    :param a2: coefficient
    :param a3: coefficient
    """

    a0: float
    a1: float
    a2: float
    a3: float


@dataclass
class TemperatureFrequencyCoefficients:
    """
    :param g: coefficient
    :param h: coefficient
    :param i: coefficient
    :param j: coefficient
    :param f0: coefficient
    """

    g: float
    h: float
    i: float
    j: float
    f0: float


@dataclass
class PressureCoefficients:
    """
    :param pa0: coefficient
    :param pa1: coefficient
    :param pa2: coefficient
    :param ptca0: coefficient
    :param ptca1: coefficient
    :param ptca2: coefficient
    :param ptcb0: coefficient
    :param ptcb1: coefficient
    :param ptcb2: coefficient
    :param ptempa0: coefficient
    :param ptempa1: coefficient
    :param ptempa2: coefficient
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
class PressureDigiquartzCoefficients:
    """
    :param a0: coefficient
    :param a1: coefficient
    :param a2: coefficient
    :param b0: coefficient
    :param b1: coefficient
    :param b2: coefficient
    :param c0: coefficient
    :param c1: coefficient
    :param c2: coefficient
    :param d0: coefficient
    :param d1: coefficient
    :param d2: coefficient
    """

    c1: float
    c2: float
    c3: float
    d1: float
    d2: float
    t1: float
    t2: float
    t3: float
    t4: float
    t5: float
    AD590M: float
    AD590B: float

    def __init__(
        self,
        c1: float,
        c2: float,
        c3: float,
        d1: float,
        d2: float,
        t1: float,
        t2: float,
        t3: float,
        t4: float,
        t5: float,
        AD590M: Optional[float] = None,
        AD590B: Optional[float] = None,
    ):
        self.c1 = c1
        self.c2 = c2
        self.c3 = c3
        self.d1 = d1
        self.d2 = d2
        self.t1 = t1
        self.t2 = t2
        self.t3 = t3
        self.t4 = t4
        self.t5 = t5
        if AD590M is not None:
            self.AD590M = AD590M
        else:
            self.AD590M = 0.028927  # default value for SBE16plus V2 and SBE19plus V2

        if AD590B is not None:
            self.AD590B = AD590B
        else:
            self.AD590B = -41.1733  # default value for SBE16plus V2 and SBE19plus V2


@dataclass
class ConductivityCoefficients:
    """
    :param g: coefficient
    :param h: coefficient
    :param i: coefficient
    :param j: coefficient
    :param cpcor: compressibility coefficient
    :param ctcor: coefficient
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
    :param a0: coefficient
    :param a1: coefficient
    :param a2: coefficient
    :param b0: coefficient
    :param b1: coefficient
    :param c0: coefficient
    :param c1: coefficient
    :param c2: coefficient
    :param e: coefficient
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
    :param ta0: coefficient
    :param ta1: coefficient
    :param ta2: coefficient
    :param ta3: coefficient
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
    :param a: coefficient
    :param b: coefficient
    :param c: coefficient
    :param e: coefficient
    :param d0: tau(T,P) coefficient
    :param d1: tau(T,P) coefficient
    :param d2: tau(T,P) coefficient
    :param h1: hysteresis correction coefficient
    :param h2: hysteresis correction coefficient
    :param h3: hysteresis correction coefficient
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
    :param multiplier: 1.0 for units of μmol photons/m2*s
    """

    im: float
    a0: float
    a1: float
    multiplier: float


@dataclass
class SPARCoefficients:
    """
    :param im: immersion coefficient
    :param a0: calibration slope
    :param a1: calibration offset
    :param conversion_factor: 1.0 for units of μmol photons/m2*s
    """

    im: float
    a0: float
    a1: float
    conversion_factor: float


@dataclass
class TemperatureSeaFETCoefficients:
    """
    :param tdfa0: temeperature coefficient
    :param tdfa1: temeperature coefficient
    :param tdfa2: temeperature coefficient
    :param tdfa3: temeperature coefficient
    """

    tdfa0: float
    tdfa1: float
    tdfa2: float
    tdfa3: float


class PHSeaFETInternalCoefficients:
    """
    :param kdf0: internal K0 coefficient
    :param kdf2: internal K2 coefficient

    :param k0: Deprecated, use kdf0
    :param k2: Deprecated, use kdf2
    :param int_k0: Deprecated, use kdf0
    :param int_k2: Deprecated, use kdf2
    """

    def __init__(
        self,
        kdf0: float = 0,
        kdf2: float = 0,
        k0: Optional[float] = None,
        k2: Optional[float] = None,
        int_k0: Optional[float] = None,
        int_k2: Optional[float] = None,
    ):
        self.kdf0 = kdf0
        self.kdf2 = kdf2
        if k0 is not None:
            self.k0 = k0
            self.kdf0 = k0
            warnings.warn("k0 is deprecated, use kdf0", DeprecationWarning)
        if k2 is not None:
            self.k2 = k2
            self.kdf2 = k2
            warnings.warn("k2 is deprecated, use kdf2", DeprecationWarning)
        if int_k0 is not None:
            self.int_k0 = int_k0
            self.kdf0 = int_k0
            warnings.warn("int_k0 is deprecated, use kdf0", DeprecationWarning)
        if int_k2 is not None:
            self.int_k2 = int_k2
            self.kdf2 = int_k2
            warnings.warn("int_k2 is deprecated, use kdf2", DeprecationWarning)


class PHSeaFETExternalCoefficients:
    """
    :param k0: external K0 coefficient
    :param k2: external K2 coefficient
    :param k2_poly_order: order of K2 pressure compensation calculation
    :param k2f0: f(P) coefficient
    :param k2f1: f(P) coefficient
    :param k2f2: f(P) coefficient
    :param k2f3: f(P) coefficient
    :param fp_poly_order: order of pressure compensation calculation
    :param f0: f(P) coefficient
    :param f1: f(P) coefficient
    :param f2: f(P) coefficient
    :param f3: f(P) coefficient
    :param f4: f(P) coefficient
    :param f5: f(P) coefficient
    :param f6: f(P) coefficient
    """

    def __init__(
        self,
        k0: float = 0,
        k2: float = 0,
        k2_poly_order=0,
        k2f0: float = 0,
        k2f1: float = 0,
        k2f2: float = 0,
        k2f3: float = 0,
        fp_poly_order: float = 0,
        f0: float = 0,
        f1: float = 0,
        f2: float = 0,
        f3: float = 0,
        f4: float = 0,
        f5: float = 0,
        f6: float = 0,
        ext_k0: Optional[float] = None,
        ext_k2: Optional[float] = None,
    ):
        self.k0 = k0
        self.k2 = k2
        self.k2_poly_order = k2_poly_order
        self.k2f0 = k2f0
        self.k2f1 = k2f1
        self.k2f2 = k2f2
        self.k2f3 = k2f3
        self.fp_poly_order = fp_poly_order
        self.f0 = f0
        self.f1 = f1
        self.f2 = f2
        self.f3 = f3
        self.f4 = f4
        self.f5 = f5
        self.f6 = f6

        if ext_k0 is not None:
            self.ext_k0 = ext_k0
            self.k0 = ext_k0
            warnings.warn("ext_k0 is deprecated, use k0", DeprecationWarning)
        if ext_k2 is not None:
            self.ext_k2 = ext_k2
            self.k2 = ext_k2
            warnings.warn("ext_k2 is deprecated, use k2", DeprecationWarning)


@dataclass
class AltimeterCoefficients:
    """
    :param slope: coefficient
    :param offset: coefficient
    """

    slope: float
    offset: float
