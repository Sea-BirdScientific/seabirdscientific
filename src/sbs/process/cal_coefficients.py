#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""Some cal data for individual serial number units.

TODO: move outside of python code?"""

# Native imports
from dataclasses import dataclass

# Third-party imports

# Sea-Bird imports

# Internal imports


@dataclass
class SN6130: # cal coefficients for E8001
    # temperature
    A0 = 1.28015621e-003
    A1 = 2.58367774e-004
    A2 = -1.39527596e-008
    A3 = 1.39024630e-007

    # pressure 
    PA0 = 6.51546669e-002
    PA1 = 1.54424116e-003
    PA2 = 6.13653149e-012
    PTCA0 = 5.24108391e+005
    PTCA1 = 5.47371611e+000
    PTCA2 = -1.53365246e-001
    PTCB0 = 2.59008750e+001
    PTCB1 = 7.75000000e-004
    PTCB2 = 0.00000000e+000
    PTEMPA0 = -6.28624239e+001
    PTEMPA1 = 5.41620689e+001
    PTEMPA2 = -2.96026659e-001

    # conductivity
    G = -1.02365730e+00
    H = 1.49306223e-01
    I = -5.63206187e-04
    J = 6.51821272e-05
    CPCOR = -9.570000e-08
    CTCOR = 3.250000e-06
    WBOTC = 0.00000000e+000


@dataclass
class SN03706385: # cal coefficients for SBE37SM-6385
    # temperature
    A0 = 1.62680600e-004
    A1 = 2.45695900e-004
    A2 = -2.34056800e-007
    A3 = 9.32450300e-008

    # pressure 
    PA0 = 6.60137700e-001
    PA1 = 9.63297700e-003
    PA2 = -1.65805500e-010
    PTCA0 = 5.24536600e+005
    PTCA1 = 4.77095200e+000
    PTCA2 = -9.34045100e-002
    PTCB0 = 2.60416300e+001
    PTCB1 = -1.67500000e-003
    PTCB2 = 0.00000000e+000
    PTEMPA0 = -6.83823200e+001
    PTEMPA1 = 5.19557000e-002
    PTEMPA2 = -6.69832900e-007

    # conductivity
    G = -9.99850800e-001
    H = 1.26606700e-001
    I = -2.89934500e-004
    J = 3.54390800e-005
    CPCOR = -9.570000e-08
    CTCOR = 3.250000e-06
    WBOTC = 1.47950000e-006


@dataclass
class SN03716125:
    # temperature
    A0 = -1.89840900e-004
    A1 = 3.13608700e-004
    A2 = -4.46184800e-006
    A3 = 2.01978100e-007

    # conductivity
    G = -1.00193200e+000
    H = 1.28425000e-001
    I = -9.96045200e-005
    J = 2.34704100e-005
    CPCOR = -9.57000000e-008
    CTCOR = 3.2500e-006
    WBOTC = -1.49105000e-007

    # pressure
    PA0 = 1.05494900e-003
    PA1 = 4.85370000e-004
    PA2 = -2.07071300e-012
    PTEMPA0 = -6.15215600e+001
    PTEMPA1 = 5.18166200e-002
    PTEMPA2 = -2.89409800e-007
    PTCA0 = 5.23607500e+005
    PTCA1 = 4.34649000e+000
    PTCA2 = -1.75847200e-001
    PTCB0 = 2.52622500e+001
    PTCB1 = 6.50000000e-004
    PTCB2 = 0.00000000e+000

@dataclass
class SN06302568:
    # Oxygen
    a0 = 1.0513
    a1 = -1.5e-3
    a2 = 4.1907e-1
    b0 = -2.5004e-1
    b1 = 1.6524
    c0 = 1.0355e-1
    c1 = 4.4295e-3
    c2 = 6.0011e-5
    e = 1.1e-2