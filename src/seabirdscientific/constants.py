"""Physical constants and unit conversion factors shared across the toolkit."""

#: Counts per volt for voltage channels in hex data (counts / COUNTS_TO_VOLTS = volts)
COUNTS_TO_VOLTS = 13107
#: psi per decibar (dbar * DBAR_TO_PSI = psi)
DBAR_TO_PSI = 1.450377
#: [Coulombs mol^{-1}] Faraday constant from SBS application note 99
F = 96485.365
#: Value that marks bad or missing samples, which are flagged rather than removed
FLAG_VALUE = -9.99e-29
#: Meters of fresh water per decibar, for depth from pressure in fresh water
FRESHWATER_PRESSURE_TO_DEPTH = 1.019716
#: IPTS-68 temperature per ITS-90 temperature (T68 = T90 * ITS90_TO_IPTS68), taken from
#: https://blog.seabird.com/ufaqs/what-is-the-difference-in-temperature-expressions-between-ipts-68-and-its-90/
ITS90_TO_IPTS68 = 1.00024
#: [K] 0 degrees C in kelvin
KELVIN_OFFSET_0C = 273.15
#: [K] 25 degrees C in kelvin
KELVIN_OFFSET_25C = 298.15
#: Feet per meter
METERS_TO_FEET = 3.28084
#: Dissolved oxygen in mg/L per mL/L
OXYGEN_MLPERL_TO_MGPERL = 1.42903
#: Dissolved oxygen in umol/m^3 per mL/L (44.66 umol/L x 1000 L/m^3); divide by seawater density
#: in kg/m^3 to get umol/kg
OXYGEN_MLPERL_TO_UMOLPERKG = 44660
#: Dissolved oxygen in umol/L per mL/L
OXYGEN_MLPERL_TO_UMOLPERL = 44.66
#: SBE 63 raw oxygen phase per volt (phase / OXYGEN_PHASE_TO_VOLTS = volts), from the manual
OXYGEN_PHASE_TO_VOLTS = 39.457071
#: Decibars per psi (psi * PSI_TO_DBAR = dbar)
PSI_TO_DBAR = 0.6894759
#: [J K^{-1} mol^{-1}] Gas constant from SBS application note 99
R = 8.3144621
#: [dbar] Atmospheric pressure at sea level
SEA_LEVEL_PRESSURE_DBAR = 10.1325
#: [psi] Atmospheric pressure at sea level, for converting between absolute and gauge pressure
SEA_LEVEL_PRESSURE_PSI = 14.7
#: Seconds from the Unix epoch (1970-01-01) to 2000-01-01, the reference for instrument timestamps
SECONDS_BETWEEN_EPOCH_AND_2000 = 946684800
#: Temperature coefficient of conductivity (fraction per degree C), used to normalize specific
#: conductance to 25 degrees C
THERMAL_CONDUCTIVITY_COEFF = 0.02
#: micro moles of nitrate to milligrams of nitrogen per liter
UMNO3_TO_MGNL = 0.014007
