from seabirdscientific.cal_coefficients import (
    ConductivityCoefficients,
    Oxygen43Coefficients,
    Oxygen63Coefficients,
    PARCoefficients,
    PH18Coefficients,
    PHSeaFETInternalCoefficients,
    PHSeaFETExternalCoefficients,
    PressureCoefficients,
    PressureDigiquartzCoefficients,
    TemperatureCoefficients,
    TemperatureFrequencyCoefficients,
    Thermistor63Coefficients,
    ECOCoefficients,
)


# @dataclass
# class SN6130:  # cal coefficients for example 19pV2 data
# temperature
temperature_coefs_sn6130 = TemperatureCoefficients(
    a0=1.28015621e-003,
    a1=2.58367774e-004,
    a2=-1.39527596e-008,
    a3=1.39024630e-007,
)

temperature_coefs_sn5102 = TemperatureFrequencyCoefficients(
    g=4.39377604e-003, h=6.43356818e-004, i=2.23495143e-005, j=2.02840480e-006, f0=1000.000
)

temperature_coefs_sn5109 = TemperatureFrequencyCoefficients(
    g=4.39126107e-003, h=6.42119447e-004, i=2.28636318e-005, j=2.17029078e-006, f0=1000.000
)

temperature_coefs_sn03906502 = TemperatureCoefficients(
    a0=-0.0008170231,
    a1=0.0003358073,
    a2=-0.000006115756,
    a3=2.096198e-7,
)

pressure_coefs_sn03906502 = PressureCoefficients(
    pa0=0.5459716,
    pa1=0.0009598281,
    pa2=-2.781709e-12,
    ptca0=2018.425,
    ptca1=87.51056,
    ptca2=-2.633389,
    ptcb0=25.02572,
    ptcb1=0.0003508772,
    ptcb2=0,
    ptempa0=685.9287,
    ptempa1=0.005264199,
    ptempa2=-0.000001318757,
)

pressure_coefs_sn6130 = PressureCoefficients(
    pa0=6.51546669e-002,
    pa1=1.54424116e-003,
    pa2=6.13653149e-012,
    ptca0=5.24108391e005,
    ptca1=5.47371611e000,
    ptca2=-1.53365246e-001,
    ptcb0=2.59008750e001,
    ptcb1=7.75000000e-004,
    ptcb2=0.00000000e000,
    ptempa0=-6.28624239e001,
    ptempa1=5.41620689e001,
    ptempa2=-2.96026659e-001,
)

pressure_coefs_sn0936 = PressureDigiquartzCoefficients(
    c1=-4.164639e004,
    c2=-5.769818e-001,
    c3=1.259640e-002,
    d1=3.483300e-002,
    d2=0.000000e000,
    t1=3.004422e001,
    t2=-4.702082e-004,
    t3=4.039850e-006,
    t4=3.117530e-009,
    t5=0,
    AD590M=1.281400e-002,
    AD590B=-9.348340e000,
)

conductivity_coefs_sn6130 = ConductivityCoefficients(
    g=-1.02365730e00,
    h=1.49306223e-01,
    i=-5.63206187e-04,
    j=6.51821272e-05,
    cpcor=-9.570000e-08,
    ctcor=3.250000e-06,
    wbotc=0.00000000e000,
)

conductivity_coefs_sn3569 = ConductivityCoefficients(
    g=-9.99280446e000,
    h=1.20387434e000,
    i=-1.12890185e-003,
    j=1.30029613e-004,
    cpcor=-9.57000000e-008,
    ctcor=3.2500e-006,
    wbotc=0,
)

conductivity_coefs_sn2206 = ConductivityCoefficients(
    g=-1.00554193e001,
    h=1.36309522e000,
    i=1.92287707e-004,
    j=5.92286921e-005,
    cpcor=-9.57000000e-008,
    ctcor=3.2500e-006,
    wbotc=0,
)

# ECO
chlorophyll_a_coefs_sn6130 = ECOCoefficients(
    slope=10,
    offset=0.0680,
)

chlorophyll_a_coefs_sn3122 = ECOCoefficients(
    slope=13,
    offset=0.0380,
)


# cal coefficients for SBE37SM-6385
temperature_coefs_sn6385 = TemperatureCoefficients(
    a0=1.62680600e-004,
    a1=2.45695900e-004,
    a2=-2.34056800e-007,
    a3=9.32450300e-008,
)

pressure_coefs_sn6385 = PressureCoefficients(
    pa0=6.60137700e-001,
    pa1=9.63297700e-003,
    pa2=-1.65805500e-010,
    ptca0=5.24536600e005,
    ptca1=4.77095200e000,
    ptca2=-9.34045100e-002,
    ptcb0=2.60416300e001,
    ptcb1=-1.67500000e-003,
    ptcb2=0.00000000e000,
    ptempa0=-6.83823200e001,
    ptempa1=5.19557000e-002,
    ptempa2=-6.69832900e-007,
)

conductivity_coefs_sn6385 = ConductivityCoefficients(
    g=-9.99850800e-001,
    h=1.26606700e-001,
    i=-2.89934500e-004,
    j=3.54390800e-005,
    cpcor=-9.570000e-08,
    ctcor=3.250000e-06,
    wbotc=1.47950000e-006,
)


# @dataclass
# cal coefficients for SN03716125:
temperature_coefs_sn16125 = TemperatureCoefficients(
    a0=-1.89840900e-004,
    a1=3.13608700e-004,
    a2=-4.46184800e-006,
    a3=2.01978100e-007,
)

conductivity_coefs_sn16125 = ConductivityCoefficients(
    g=-1.00193200e000,
    h=1.28425000e-001,
    i=-9.96045200e-005,
    j=2.34704100e-005,
    cpcor=-9.57000000e-008,
    ctcor=3.2500e-006,
    wbotc=-1.49105000e-007,
)

pressure_coefs_sn16125 = PressureCoefficients(
    pa0=1.05494900e-003,
    pa1=4.85370000e-004,
    pa2=-2.07071300e-012,
    ptca0=5.23607500e005,
    ptca1=4.34649000e000,
    ptca2=-1.75847200e-001,
    ptcb0=2.52622500e001,
    ptcb1=6.50000000e-004,
    ptcb2=0.00000000e000,
    ptempa0=-6.15215600e001,
    ptempa1=5.18166200e-002,
    ptempa2=-2.89409800e-007,
)


# cal coefficients for SN06302568
oxygen_63_coefs_sn2568 = Oxygen63Coefficients(
    a0=1.0513,
    a1=-1.5e-3,
    a2=4.1907e-1,
    b0=-2.5004e-1,
    b1=1.6524,
    c0=1.0355e-1,
    c1=4.4295e-3,
    c2=6.0011e-5,
    e=1.1e-2,
)

thermistor_63_coefs_sn2568 = Thermistor63Coefficients(
    ta0=7.059180e-4,
    ta1=2.504670e-4,
    ta2=7.402389e-7,
    ta3=9.756123e-8,
)


# cal coefficients for SN3287
oxygen_43_coefs_sn3287 = Oxygen43Coefficients(
    soc=0.4792,
    v_offset=-0.484,
    tau_20=2.05,
    a=-3.661e-3,
    b=1.7812e-4,
    c=-2.5951e-6,
    e=0.036,
    d0=0,
    d1=1.92634e-4,
    d2=-4.64803e-2,
    h1=-3.3e-2,
    h2=5e3,
    h3=1.45e3,
)


# cal coefficients for SN431686
oxygen_43_coefs_sn1686 = Oxygen43Coefficients(
    soc=3.898e-1,
    v_offset=-0.498,
    tau_20=1.08,
    a=-3.3982e-3,
    b=1.5817e-4,
    c=-2.6651e-6,
    e=0.036,
    d0=2.5826,
    d1=1.92634e-4,
    d2=-4.64803e-2,
    h1=-3.3e-2,
    h2=5e3,
    h3=1.45e3,
)

# cal coefficients for SN1590
oxygen_43_coefs_sn1590 = Oxygen43Coefficients(
    soc=5.7640e-001,
    v_offset=-0.5278,
    tau_20=1.5700,
    a=-2.1585e-003,
    b=1.0617e-004,
    c=1.9903e-006,
    e=3.6000e-002,
    d0=2.5826e000,
    d1=1.92634e-004,
    d2=-4.64803e-002,
    h1=-3.3000e-002,
    h2=5.0000e003,
    h3=1.4500e003,
)

# cal coefficients for SN0680
oxygen_43_coefs_sn0680 = Oxygen43Coefficients(
    soc=4.9197e-001,
    v_offset=-0.4928,
    tau_20=1.1900,
    a=-2.7974e-003,
    b=1.2124e-004,
    c=-1.9448e-006,
    e=3.6000e-002,
    d0=2.5826e000,
    d1=1.92634e-004,
    d2=-4.64803e-002,
    h1=-3.3000e-002,
    h2=5.0000e003,
    h3=1.4500e003,
)

# cal coefficients for SN180726
ph_coefs_sn0762 = PH18Coefficients(
    slope=4.6067,
    offset=2.5157,
)

# cal coefficients for SN0709
ph_coefs_sn0709 = PH18Coefficients(
    slope=4.4559,
    offset=2.5700,
)

# cal coefficients for SN0411
par_coefs_sn0411 = PARCoefficients(
    im=1.3589,
    a0=1.372,
    a1=0.8839,
    multiplier=1,
)

ph_seafet_internal_coeffs = PHSeaFETInternalCoefficients(
    k0=-1.354051,
    k2=-0.001272858,
)

ph_seafet_external_coeffs = PHSeaFETExternalCoefficients(
    k0=-1.345407,
    k2=-0.001106957,
    f1=0,
    f2=0,
    f3=0,
    f4=0,
    f5=0,
    f6=0,
)
