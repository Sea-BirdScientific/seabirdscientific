#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""TODO: NEW interpret_sbs_variable docstring"""

"""
DESCRIPTION:
Converts the strings that define SBE variables to strings that matlab can
use as functions. Gives units for each variable and defines whether
it is a number or a text string. 

All variables imported from SBEDataProcessing_7.26.4.pdf. There may be
additional variables added in future versions of Sea-Bird Data processing
modules that are not included here. Please email kmartini@seabird.com if
you find one to keep the code current.

INPUT:
  sbs_var         =   text string with variable name

OUTPUT: 
  kvar_name       =   converted, matlab code readable variable name
  kvar_format     =   defines whether variable is a float or a string
  kvar_units      =   variable units. Empty string if undefined.


KiM MARTiNi 07.2017
Sea-Bird Scientific 
"""

# Native imports
import re

# Third-party imports

# Sea-Bird imports

# Internal imports


def interpret_sbs_variable(sbs_var):
    """Converts the strings that define SBE variables to strings that matlab can
    use as functions. Gives units for each variable and defines whether
    it is a number or a text string.

    Args:
        sbs_var (string): Name of the dataset

    Returns:
        list[str]: The name with special characters removed, data format, and the units.
    """
    # TEMPERATURE
    # primary sensor
    kvar_format = r"%f"

    if sbs_var in {"t090Cm", "t4990C", "tnc90C", "tv290C", "t090C"}:
        kvar_name = "t090C"
        kvar_units = "ITS-90, deg C"
    elif sbs_var in {"t090F", "t4990F", "tnc90F", "tv290F"}:
        kvar_name = "t090F"
        kvar_units = "ITS-90, deg F"
    elif sbs_var in {"t068C", "t4968C", "tnc68C", "tv268C"}:
        kvar_name = "t068C"
        kvar_units = "ITS-68, deg C"
    elif sbs_var in {"t068F", "t4968F", "tnc68F", "tv268F"}:
        kvar_name = "t068F"
        kvar_units = "ITS-68, deg F"
    elif sbs_var in {"t090"}:
        kvar_name = "t090"
        kvar_units = "SBE49, deg C"

    # secondary temperature
    elif sbs_var in {"t190C", "tnc290C"}:
        kvar_name = "t190C"
        kvar_units = "ITS-90, deg C"
    elif sbs_var in {"t190F", "tnc290F"}:
        kvar_name = "t190F"
        kvar_units = "ITS-90, deg F"
    elif sbs_var in {"t168C", "tnc268C"}:
        kvar_name = "t168C"
        kvar_units = "ITS-68, deg C"
    elif sbs_var in {"t168F", "tnc268F"}:
        kvar_name = "t168F"
        kvar_units = "ITS-68, deg F"

    # SBE38 REFERENCE TEMPERATURE
    # primary sensor
    elif sbs_var in {"t3890C", "t38_90C"}:
        kvar_name = "t03890C"
        kvar_units = "ITS-90, deg C"
    elif sbs_var in {"t3890F", "t38_90F"}:
        kvar_name = "t03890F"
        kvar_units = "ITS-90, deg F"
    elif sbs_var in {"t3868C", "t38_68C"}:
        kvar_name = "t03868C"
        kvar_units = "ITS-68, deg C"
    elif sbs_var in {"t3868F", "t38_68F"}:
        kvar_name = "t03868F"
        kvar_units = "ITS-68, deg F"
    # secondary sensor
    elif sbs_var in {"t3890C1"}:
        kvar_name = "t13890C"
        kvar_units = "ITS-90, deg C"
    elif sbs_var in {"t3890F1"}:
        kvar_name = "t13890F"
        kvar_units = "ITS-90, deg F"
    elif sbs_var in {"t3868C1"}:
        kvar_name = "t13868C"
        kvar_units = "ITS-68, deg C"
    elif sbs_var in {"t3868F1"}:
        kvar_name = "t13868F"
        kvar_units = "ITS-68, deg F"

    # TEMPERATURE DIFFERENCE
    elif sbs_var in {"T2-T190C"}:
        kvar_name = "tdiff90C"
        kvar_units = "ITS-90, deg C"
    elif sbs_var in {"T2-T190F"}:
        kvar_name = "tdiff90F"
        kvar_units = "ITS-90, deg F"
    elif sbs_var in {"T2-T168C"}:
        kvar_name = "tdiff68C"
        kvar_units = "IPTS-68, deg C"
    elif sbs_var in {"T2-T168F"}:
        kvar_name = "tdiff68F"
        kvar_units = "IPTS-68, deg F"

    # CONDUCTIVITY
    # primary conductivity
    elif sbs_var in {"c_S/m", "cond0S/m", "c0S/m"}:
        kvar_name = "c0Sm"
        kvar_units = "S/m"
    elif sbs_var in {"c_mS/cm", "cond0mS/cm", "c0mS/cm"}:
        kvar_name = "c0mScm"
        kvar_units = "mS/cm"
    elif sbs_var in {"c_uS/cm", "cond0uS/cm", "c0uS/cm"}:
        kvar_name = "c0uScm"
        kvar_units = "microS/cm"
    # secondary conductivity
    elif sbs_var in {"c1S/m"}:
        kvar_name = "c1Sm"
        kvar_units = "S/m"
    elif sbs_var in {"c1mS/cm"}:
        kvar_name = "c1mScm"
        kvar_units = "mS/cm"
    elif sbs_var in {"c1uS/cm"}:
        kvar_name = "c1uScm"
        kvar_units = "microS/cm"

    # CONDUCTIVITY DIFFERENCE
    elif sbs_var in {"C2-C1S/m"}:
        kvar_name = "cdiffSm"
        kvar_units = "S/m"
    elif sbs_var in {"C2-C1mS/cm"}:
        kvar_name = "cdiffmScm"
        kvar_units = "mS/cm"
    elif sbs_var in {"C2-C1uS/cm"}:
        kvar_name = "cdiffuScm"
        kvar_units = "microS/cm"

    # PRESSURE
    # user entered pressure values
    elif sbs_var in {"pr"}:  # SBE 49
        kvar_name = "p"
        kvar_units = "db"
    elif sbs_var in {"prM"}:
        kvar_name = "pM"
        kvar_units = "db (user entered)"
    elif sbs_var in {"prE"}:
        kvar_name = "pE"
        kvar_units = "psi (user entered)"
    # digiquartz pressure
    elif sbs_var in {"prDM"}:
        kvar_name = "pm"
        kvar_units = "db"
    elif sbs_var in {"prDE"}:
        kvar_name = "pe"
        kvar_units = "psi"
    # strain gauge pressure
    elif sbs_var in {"prSM", "prdM"}:
        kvar_name = "pm"
        kvar_units = "db"
    elif sbs_var in {"prSE", "prdE"}:
        kvar_name = "pe"
        kvar_units = "psi"
    # pressure temperature
    elif sbs_var in {"ptempC"}:
        kvar_name = "ptempC"
        kvar_units = "deg C"
    elif sbs_var in {"ptempF"}:
        kvar_name = "ptempF"
        kvar_units = "deg F"
    # FGP Pressure
    elif sbs_var in {"fgp0"}:
        kvar_name = "fgp0"
        kvar_units = "KPa"
    elif sbs_var in {"fgp1"}:
        kvar_name = "fgp1"
        kvar_units = "KPa"
    elif sbs_var in {"fgp2"}:
        kvar_name = "fgp2"
        kvar_units = "KPa"
    elif sbs_var in {"fgp3"}:
        kvar_name = "fgp3"
        kvar_units = "KPa"
    elif sbs_var in {"fgp4"}:
        kvar_name = "fgp4"
        kvar_units = "KPa"
    elif sbs_var in {"fgp5"}:
        kvar_name = "fgp5"
        kvar_units = "KPa"
    elif sbs_var in {"fgp6"}:
        kvar_name = "fgp6"
        kvar_units = "KPa"
    elif sbs_var in {"fgp7"}:
        kvar_name = "fgp7"
        kvar_units = "KPa"
    # SBE 50
    elif sbs_var in {"pr50M"}:
        kvar_name = "p050m"
        kvar_units = "db"
    elif sbs_var in {"pr50E"}:
        kvar_name = "p050e"
        kvar_units = "psi"
    elif sbs_var in {"pr50M1"}:
        kvar_name = "p150m"
        kvar_units = "db"
    elif sbs_var in {"pr50E1"}:
        kvar_name = "p150e"
        kvar_units = "psi"

    # FREQUENCY OUTPUT
    elif sbs_var in {
        "f0",
        "f1",
        "f2",
        "f3",
        "f4",
        "f5",
        "f6",
        "f7",
        "f8",
        "f9",
        "f10",
        "f11",
        "f12",
        "f13",
        "f14",
        "f15",
        "f16",
        "f17",
        "f18",
        "f19",
        "f20",
        "f21",
        "f22",
        "f23",
        "f24",
        "f25",
        "f26",
        "f27",
        "f28",
        "f29",
        "f30",
        "f31",
        "f32",
        "f33",
        "f34",
        "f35",
        "f36",
    }:
        kvar_name = sbs_var
        kvar_units = "Hz"

    # VOLTAGE OUTPUT
    elif sbs_var in {
        "v0",
        "v1",
        "v2",
        "v3",
        "v4",
        "v5",
        "v6",
        "v7",
        "v8",
        "v9",
        "v10",
        "v11",
        "v12",
        "v13",
        "v14",
        "v15",
    }:
        kvar_name = sbs_var
        kvar_units = "V"

    # GTD-DO SENSORS
    # PRESSURE
    elif sbs_var in {"GTDDOP0", "GTDDOP1", "GTDDOPdiff"}:
        kvar_name = sbs_var
        kvar_units = "mb"
    # TEMPERATURE
    elif sbs_var in {"GTDDOT0", "GTDDOT1", "GTDDOTdiff"}:
        kvar_name = sbs_var
        kvar_units = "deg C"
    # GTD-N2 SENSOR
    # PRESSURE
    elif sbs_var in {"GTDN2P0", "GTDN2P1", "GTDN2Pdiff"}:
        kvar_name = sbs_var
        kvar_units = "mb"
    # TEMPERATURE
    elif sbs_var in {"GTDN2T0", "GTDN2T1", "GTDN2Tdiff"}:
        kvar_name = sbs_var
        kvar_units = "deg C"

    # TIME
    # elapsed time
    elif sbs_var in {"timeS"}:
        kvar_name = "timeS"
        kvar_units = "elapsed time seconds"
    elif sbs_var in {"timeM"}:
        kvar_name = "timeM"
        kvar_units = "elapsed time minutes"
    elif sbs_var in {"timeH"}:
        kvar_name = "timeH"
        kvar_units = "elapsed time hours"
    elif sbs_var in {"timeJ"}:
        kvar_name = "timeJ"
        kvar_units = "elapsed time Julian Days"
    elif sbs_var in {"timeN"}:
        kvar_name = "timeN"
        kvar_units = "From NMEA: seconds since January 1, 2000"
    elif sbs_var in {"timeQ"}:
        kvar_name = "timeQ"
        kvar_units = "From NMEA: seconds since January 1, 1970"
    elif sbs_var in {"timeK"}:
        kvar_name = "timeK"
        kvar_units = "SBE timestamp: seconds since January 1, 2000"
    elif sbs_var in {"timeJV2"}:  # julian days
        kvar_name = "jdays"
        kvar_units = "SBE timestamp: julian days"
    elif sbs_var in {"timeSCP"}:  # julian days
        kvar_name = "jdays"
        kvar_units = "SBE timestamp: julian days"
    elif sbs_var in {"timeY"}:
        kvar_name = "timeY"
        kvar_units = "Computer Time: seconds since January 1, 1970"

    # POSITION
    elif sbs_var in {"latitude"}:
        kvar_name = "lat"
        kvar_units = "from NMEA degrees"
    elif sbs_var in {"longitude"}:
        kvar_name = "lon"
        kvar_units = "from NMEA degrees"

    # DEPTH
    elif sbs_var in {"depSM"}:
        kvar_name = "depSM"
        kvar_units = "depth salt water, m"
    elif sbs_var in {"depSF"}:
        kvar_name = "depSF"
        kvar_units = "depth salt water, ft"
    elif sbs_var in {"depFM"}:
        kvar_name = "depFM"
        kvar_units = "depth fresh water, m"
    elif sbs_var in {"depFF"}:
        kvar_name = "depSM"
        kvar_units = "depth fresh water, ft"
    elif sbs_var in {"dNMEA"}:
        kvar_name = "depSM"
        kvar_units = "depth from NMEA salt water, m"

    # POTENTIAL TEMPERATURE
    # I use the abbreviation of the greek symbol theta to designate
    # potential temperature
    elif sbs_var in {"potemp090C"}:
        kvar_name = "th090C"
        kvar_units = "ITS-90, deg C"
    elif sbs_var in {"potemp090F"}:
        kvar_name = "th090F"
        kvar_units = "ITS-90, deg F"
    elif sbs_var in {"potemp068C"}:
        kvar_name = "th068C"
        kvar_units = "ITS-68, deg C"
    elif sbs_var in {"potemp068F"}:
        kvar_name = "th068F"
        kvar_units = "ITS-68, deg F"
    elif sbs_var in {"potemp190C"}:
        kvar_name = "th190C"
        kvar_units = "ITS-90, deg C"
    elif sbs_var in {"potemp190F"}:
        kvar_name = "th190F"
        kvar_units = "ITS-90, deg F"
    elif sbs_var in {"potemp168C"}:
        kvar_name = "th168C"
        kvar_units = "ITS-68, deg C"
    elif sbs_var in {"potemp168F"}:
        kvar_name = "th168F"
        kvar_units = "ITS-68, deg F"

    # POTENTIAL TEMPEATURE DIFFERENCE
    elif sbs_var in {"potemp90Cdiff"}:
        kvar_name = "th090Cdiff"
        kvar_units = "ITS-90, deg C"
    elif sbs_var in {"potemp90Fdiff"}:
        kvar_name = "th090Fdiff"
        kvar_units = "ITS-90, deg F"
    elif sbs_var in {"potemp68Cdiff"}:
        kvar_name = "th068Cdiff"
        kvar_units = "ITS-68, deg C"
    elif sbs_var in {"potemp68Fdiff"}:
        kvar_name = "th068Fdiff"
        kvar_units = "ITS-68, deg F"

    # POTENTIAL TEMPEATURE ANOMALY
    elif sbs_var in {"tha090C"}:
        kvar_name = "tha090C"
        kvar_units = "ITS-90, deg C"
    elif sbs_var in {"tha090F"}:
        kvar_name = "tha090F"
        kvar_units = "ITS-90, deg F"
    elif sbs_var in {"tha068C"}:
        kvar_name = "tha068C"
        kvar_units = "ITS-68, deg C"
    elif sbs_var in {"tha068F"}:
        kvar_name = "tha068F"
        kvar_units = "ITS-68, deg F"
    elif sbs_var in {"tha190C"}:
        kvar_name = "tha190C"
        kvar_units = "ITS-90, deg C"
    elif sbs_var in {"tha190F"}:
        kvar_name = "tha190F"
        kvar_units = "ITS-90, deg F"
    elif sbs_var in {"tha168C"}:
        kvar_name = "tha168C"
        kvar_units = "ITS-68, deg C"
    elif sbs_var in {"tha168F"}:
        kvar_name = "tha168F"
        kvar_units = "ITS-68, deg F"
    elif sbs_var in {"pta090C"}:
        kvar_name = "pta090C"
        kvar_units = "ITS-90, deg C"

    # SALINITY
    elif sbs_var in {"sal00"}:
        kvar_name = "psal0"
        kvar_units = "PSU"
    elif sbs_var in {"sal11"}:
        kvar_name = "psal1"
        kvar_units = "PSU"
    elif sbs_var in {"secS-priS"}:
        kvar_name = "psaldiff"
        kvar_units = "PSU"

    # DENSITY
    # sg: potential density, based on greek letter sg
    # sgth: potential denisty sg calculated with potential temperature
    # theta (th)
    elif sbs_var in {"density00"}:
        kvar_name = "rho0"
        kvar_units = "density, kg/m^3"
    elif sbs_var in {"sigma-�00"}:
        kvar_name = "sgth0"
        kvar_units = "sigma-theta (p=0 db), kg/m^3-1000"
    elif sbs_var in {"sigma-t00"}:
        kvar_name = "sgt0"
        kvar_units = "sigma-t (p=0 db), kg/m^3-1000"
    elif sbs_var in {"sigma-100"}:
        kvar_name = "sg10"
        kvar_units = "sigma-1 (p=1000 db), kg/m^3-1000"
    elif sbs_var in {"sigma-200"}:
        kvar_name = "sg20"
        kvar_units = "sigma-2 (p=2000 db), kg/m^3-1000"
    elif sbs_var in {"sigma-300"}:
        kvar_name = "sg30"
        kvar_units = "sigma-3 (p=3000 db), kg/m^3-1000"
    elif sbs_var in {"sigma-400"}:
        kvar_name = "sg40"
        kvar_units = "sigma-4 (p=4000 db), kg/m^3-1000"
    elif sbs_var in {"density11"}:
        kvar_name = "rho1"
        kvar_units = "density, kg/m^3"
    elif sbs_var in {"sigma-�11"}:
        kvar_name = "sgth1"
        kvar_units = "sigma-theta (p=0 db), kg/m^3-1000"
    elif sbs_var in {"sigma-t11"}:
        kvar_name = "sgt1"
        kvar_units = "sigma-t (p=0 db), kg/m^3-1000"
    elif sbs_var in {"sigma-111"}:
        kvar_name = "sg11"
        kvar_units = "sigma-1 (p=1000 db), kg/m^3-1000"
    elif sbs_var in {"sigma-211"}:
        kvar_name = "sg21"
        kvar_units = "sigma-2 (p=2000 db), kg/m^3-1000"
    elif sbs_var in {"sigma-311"}:
        kvar_name = "sg31"
        kvar_units = "sigma-3 (p=3000 db), kg/m^3-1000"
    elif sbs_var in {"sigma-411"}:
        kvar_name = "sg41"
        kvar_units = "sigma-4 (p=4000 db), kg/m^3-1000"

    # DENSITY DIFFERENCE
    elif sbs_var in {"D2-D1, d"}:
        kvar_name = "rhodiff"
        kvar_units = "density, kg/m^3"
    elif sbs_var in {"D2-D1,th"}:
        kvar_name = "sgthdiff"
        kvar_units = "sigma-theta, kg/m^3"
    elif sbs_var in {"D2-D1,t"}:
        kvar_name = "sgtdiff"
        kvar_units = "sigma-t, kg/m^3"
    elif sbs_var in {"D2-D1,1"}:
        kvar_name = "sg1diff"
        kvar_units = "sigma-1 (p=1000 db), kg/m^3"
    elif sbs_var in {"D2-D1,2"}:
        kvar_name = "sg2diff"
        kvar_units = "sigma-2 (p=1000 db), kg/m^3"
    elif sbs_var in {"D2-D1,3"}:
        kvar_name = "sg3diff"
        kvar_units = "sigma-3 (p=1000 db), kg/m^3"
    elif sbs_var in {"D2-D1,4"}:
        kvar_name = "sg4diff"
        kvar_units = "sigma-4 (p=1000 db), kg/m^3"

    # SOUND VELOCITY
    # primary sensor
    elif sbs_var in {"svCM"}:
        kvar_name = "sv0CM"
        kvar_units = "Chen-Millero, m/s"
    elif sbs_var in {"svCF"}:
        kvar_name = "sv0CF"
        kvar_units = "Chen-Millero, ft/s"
    elif sbs_var in {"svDM"}:
        kvar_name = "sv0DM"
        kvar_units = "Delgrosso, m/s"
    elif sbs_var in {"svDF"}:
        kvar_name = "sv0DF"
        kvar_units = "Delgrosso, ft/s"
    elif sbs_var in {"svWM"}:
        kvar_name = "sv0WM"
        kvar_units = "Wilson, m/s"
    elif sbs_var in {"svWF"}:
        kvar_name = "sv0WF"
        kvar_units = "Wilson, ft/s"
    # secondary sensor
    elif sbs_var in {"svCM1"}:
        kvar_name = "sv1CM"
        kvar_units = "Chen-Millero, m/s"
    elif sbs_var in {"svCF1"}:
        kvar_name = "sv1CF"
        kvar_units = "Chen-Millero, ft/s"
    elif sbs_var in {"svDM1"}:
        kvar_name = "sv1DM"
        kvar_units = "Delgrosso, m/s"
    elif sbs_var in {"svDF1"}:
        kvar_name = "sv1DF"
        kvar_units = "Delgrosso, ft/s"
    elif sbs_var in {"svWM1"}:
        kvar_name = "sv1WM"
        kvar_units = "Wilson, m/s"
    elif sbs_var in {"svWF1"}:
        kvar_name = "sv1WF"
        kvar_units = "Wilson, ft/s"
    # IOW sound velocity sensor
    elif sbs_var in {"iowSv"}:
        kvar_name = "iowSv"
        kvar_units = "IOW sound velocity sensor, m/s"
    elif sbs_var in {"sbeSv-iowSv"}:
        kvar_name = "svdiff"
        kvar_units = "SBE CTD - IOW SV sensor, m/s"
    # AVERAGE SOUND VELOCITY
    elif sbs_var in {"avgsvCM"}:
        kvar_name = "avgsvCM"
        kvar_units = "Chen-Millero, m/s"
    elif sbs_var in {"avgsvCF"}:
        kvar_name = "avgsvCF"
        kvar_units = "Chen-Millero, ft/s"
    elif sbs_var in {"avgsvDM"}:
        kvar_name = "avgsvDM"
        kvar_units = "Delgrosso, m/s"
    elif sbs_var in {"avgsvDF"}:
        kvar_name = "avgsvDF"
        kvar_units = "Delgrosso, ft/s"
    elif sbs_var in {"avgsvWM"}:
        kvar_name = "avgsvWM"
        kvar_units = "Wilson, m/s"
    elif sbs_var in {"avgsvWF"}:
        kvar_name = "avgsvWF"
        kvar_units = "Wilson, ft/s"

    # BUOYANCY
    elif sbs_var in {"N"}:
        kvar_name = "N"
        kvar_units = "cycles/hour"
    elif sbs_var in {"N^2"}:
        kvar_name = "N2"
        kvar_units = "rad^2/s^2"

    # ACCELERATION
    elif sbs_var in {"accM"}:
        kvar_name = "accM"
        kvar_units = "m/s^2"
    elif sbs_var in {"accF"}:
        kvar_name = "accF"
        kvar_units = "ft/s^2"

    # DESCENT RATE
    elif sbs_var in {"dz/dtM"}:
        kvar_name = "dzdtM"
        kvar_units = "m/s"
    elif sbs_var in {"dz/dtF"}:
        kvar_name = "dzdtF"
        kvar_units = "ft/s"

    # SEAFLOOR DEPTH
    elif sbs_var in {"sfdSM"}:
        kvar_name = "sfdSM"
        kvar_units = "salt water, m"
    elif sbs_var in {"sfdSF"}:
        kvar_name = "sfdSF"
        kvar_units = "salt water, ft"
    elif sbs_var in {"sfdFM"}:
        kvar_name = "sfdFM"
        kvar_units = "fresh water, m"
    elif sbs_var in {"sfdFF"}:
        kvar_name = "sfdFF"
        kvar_units = "fresh water, ft"

    # OTHER DERIVED VARIABLES
    # DYNAMIC METERS
    elif sbs_var in {"dm"}:
        kvar_name = "dm"
        kvar_units = "10 J/kg"
    # GEOPOTENTIAL ANOMALY
    elif sbs_var in {"gpa"}:
        kvar_name = "gpa"
        kvar_units = "J/kg"
    # PLUME ANOMALY
    elif sbs_var in {"pla"}:
        kvar_name = "pla"
        kvar_units = "?"
    # SPECIFIC VOLUME ANOMALY
    elif sbs_var in {"sva"}:
        kvar_name = "sva"
        kvar_units = "10^-8 *m^3/kg"
    # STABILITY
    elif sbs_var in {"E"}:
        kvar_name = "E"
        kvar_units = "rad^2/m"
    elif sbs_var in {"E10^-8"}:
        kvar_name = "E10e_8"
        kvar_units = "10^-8 *rad^2/m"
    # THERMOSTERIC ANOMALY
    elif sbs_var in {"tsa"}:
        kvar_name = "tsa"
        kvar_units = "10^-8 *m^3/kg"
    # SPECIFIC CONDUCTANCE
    elif sbs_var in {"specc"}:
        kvar_name = "specc"
        kvar_units = "uS/cm"
    elif sbs_var in {"speccumhoscm"}:
        kvar_name = "speccumhoscm"
        kvar_units = "umhos/cm"
    elif sbs_var in {"speccmsm"}:
        kvar_name = "speccmsm"
        kvar_units = "mS/cm"
    elif sbs_var in {"speccmmhoscm"}:
        kvar_name = "speccmmhoscm"
        kvar_units = "mmhos/cm"

    # BIOGEOCHEMICAL SENSORS

    # OXYGEN
    # RAW OXYGEN, SBE 43
    elif sbs_var in {"sbeox0V"}:
        kvar_name = "sbeox0V"
        kvar_units = "V"
    elif sbs_var in {"sbeox0F"}:
        kvar_name = "sbeox0F"
        kvar_units = "Hz"
    elif sbs_var in {"sbeox1V"}:
        kvar_name = "sbeox1V"
        kvar_units = "V"
    elif sbs_var in {"sbeox1F"}:
        kvar_name = "sbeox1F"
        kvar_units = "Hz"
    # DERIVED OXYGEN, SBE 43
    elif sbs_var in {"sbeox0ML/L"}:
        kvar_name = "sbeox0mL_L"
        kvar_units = "mL/L"
    elif sbs_var in {"sbeox0Mg/L"}:
        kvar_name = "sbeox0mg_L"
        kvar_units = "mg/L"
    elif sbs_var in {"sbeox0PS"}:
        kvar_name = "sbeox0PS"
        kvar_units = r"% saturation"
    elif sbs_var in {"sbeox0Mm/Kg"}:
        kvar_name = "sbeox0mum_kg"
        kvar_units = "micromol/kg"
    elif sbs_var in {"sbeox0Mm/L"}:
        kvar_name = "sbeox0mum_L"
        kvar_units = "micromole/L"
    elif sbs_var in {"sbeox0dOV/dT"}:
        kvar_name = "sbeox0dOV_dT"
        kvar_units = "V/s"
    elif sbs_var in {"sbeox1ML/L"}:
        kvar_name = "sbeox1mL_L"
        kvar_units = "mL/L"
    elif sbs_var in {"sbeox1Mg/L"}:
        kvar_name = "sbeox1mg_L"
        kvar_units = "mg/L"
    elif sbs_var in {"sbeox1PS"}:
        kvar_name = "sbeox1PS"
        kvar_units = r"% saturation"
    elif sbs_var in {"sbeox1Mm/Kg"}:
        kvar_name = "sbeox1mum_kg"
        kvar_units = "micromol/kg"
    elif sbs_var in {"sbeox1Mm/L"}:
        kvar_name = "sbeox1mum_L"
        kvar_units = "micromole/L"
    elif sbs_var in {"sbeox1dOV/dT"}:
        kvar_name = "sbeox1dOV_dT"
        kvar_units = "V/s"
    # Difference between SBE 43 sensors
    elif sbs_var in {"sbeox0ML/Ldiff"}:
        kvar_name = "sbeox0ML_Ldiff"
        kvar_units = "ml/l"
    elif sbs_var in {"sbeox0Mg/Ldiff"}:
        kvar_name = "sbeox0Mg_Ldiff"
        kvar_units = "mg/l"
    elif sbs_var in {"sbeox0PSdiff"}:
        kvar_name = "sbeox0PSdiff"
        kvar_units = r"% saturation"
    elif sbs_var in {"sbeox0Mm/Kgdiff"}:
        kvar_name = "sbeox0Mm_Kgdiff"
        kvar_units = "umol/kg"
    elif sbs_var in {"sbeox0Mm/Ldiff"}:
        kvar_name = "sbeox0Mm_Ldiff"
        kvar_units = "umol/l"
    # RAW OXYGEN, SBE 63
    elif sbs_var in {"sbeoxpd"}:
        kvar_name = "sbeoxpd"
        kvar_units = "usec"
    elif sbs_var in {"sbeoxpdv"}:
        kvar_name = "sbeoxpdv"
        kvar_units = "V"
    elif sbs_var in {"sbeoxpd1"}:
        kvar_name = "sbeoxpd1"
        kvar_units = "usec"
    elif sbs_var in {"sbeoxpdv1"}:
        kvar_name = "sbeoxpdv1"
        kvar_units = "sbeoxpdv1"
    elif sbs_var in {"sbeoxtv"}:
        kvar_name = "sbeoxtv"
        kvar_units = "sbeoxtv"
    elif sbs_var in {"sbeoxtv1"}:
        kvar_name = ""
        kvar_units = "sbeoxtv1"
    # RAW OXYGEN TEMPERATURE, SBE 63
    elif sbs_var in {"sbeoxTC"}:
        kvar_name = "sbeoxTC"
        kvar_units = "ITS-90, deg C"
    elif sbs_var in {"sbeoxTF"}:
        kvar_name = "sbeoxTF"
        kvar_units = "ITS-90, deg F"
    elif sbs_var in {"sbeoxTC1"}:
        kvar_name = "sbeoxTC1"
        kvar_units = "ITS-90, deg C"
    elif sbs_var in {"sbeoxTF1"}:
        kvar_name = "sbeoxTF1"
        kvar_units = "ITS-90, deg F"
    # DERIVED OXYGEN, SBE 63
    elif sbs_var in {"sbeopoxML/L"}:
        kvar_name = "sbeopoxML_L"
        kvar_units = "ml/l"
    elif sbs_var in {"sbeopoxMg/L"}:
        kvar_name = "sbeopoxMg_L"
        kvar_units = "mg/l"
    elif sbs_var in {"sbeopoxPS"}:
        kvar_name = "sbeopoxPS"
        kvar_units = r"% saturation"
    elif sbs_var in {"sbeopoxMm/Kg"}:
        kvar_name = "sbeopoxMm_Kg"
        kvar_units = "umol/kg"
    elif sbs_var in {"sbeopoxMm/L"}:
        kvar_name = "sbeopoxMm_L"
        kvar_units = "umol/l"
    elif sbs_var in {"sbeopoxML/L1"}:
        kvar_name = "sbeopoxML_L1"
        kvar_units = "ml/l"
    elif sbs_var in {"sbeopoxMg/L1"}:
        kvar_name = "sbeopoxMg_L1"
        kvar_units = "mg/l"
    elif sbs_var in {"sbeopoxPS1"}:
        kvar_name = "sbeopoxPS1"
        kvar_units = r"% saturation"
    elif sbs_var in {"SbeopoxMm/Kg1"}:
        kvar_name = "SbeopoxMm_Kg1"
        kvar_units = "umol/kg"
    elif sbs_var in {"sbeopoxMm/L1"}:
        kvar_name = "sbeopoxMm/L1"
        kvar_units = "umol/l"
    # AANDERAA OPTODE
    elif sbs_var in {"opoxML/L"}:
        kvar_name = "opoxML_L"
        kvar_units = "ml/l"
    elif sbs_var in {"opoxMg/L"}:
        kvar_name = "opoxMg_L"
        kvar_units = "mg/l"
    elif sbs_var in {"opoxPS"}:
        kvar_name = "opoxPS"
        kvar_units = r"% saturation"
    elif sbs_var in {"opoxMm/L"}:
        kvar_name = "opoxMm_L"
        kvar_units = "umol/l"
    # BECKMAN/YSI
    elif sbs_var in {"oxC"}:
        kvar_name = "oxC"
        kvar_units = "uA"
    elif sbs_var in {"oxsC"}:
        kvar_name = "oxsC"
        kvar_units = "uA"
    elif sbs_var in {"oxTC"}:
        kvar_name = "oxTC"
        kvar_units = "deg C"
    elif sbs_var in {"oxTF"}:
        kvar_name = "oxTF"
        kvar_units = "deg F"
    elif sbs_var in {"oxsTC"}:
        kvar_name = "oxsTC"
        kvar_units = "deg C"
    elif sbs_var in {"oxsTF"}:
        kvar_name = "oxsTF"
        kvar_units = "deg F"
    elif sbs_var in {"oxML/L"}:
        kvar_name = "oxML_L"
        kvar_units = "ml/l"
    elif sbs_var in {"oxMg/L"}:
        kvar_name = "oxMg_L"
        kvar_units = "mg/l"
    elif sbs_var in {"oxPS"}:
        kvar_name = "oxPS"
        kvar_units = r"% saturation"
    elif sbs_var in {"oxMm/Kg"}:
        kvar_name = "oxMm_Kg"
        kvar_units = "umol/kg"
    elif sbs_var in {"oxdOC/dT"}:
        kvar_name = "oxdOC_dT"
        kvar_units = "doc/dt"
    elif sbs_var in {"oxsML/L"}:
        kvar_name = "oxsML_L"
        kvar_units = "ml/l"
    elif sbs_var in {"oxsMg/L"}:
        kvar_name = "oxsMg_L"
        kvar_units = "mg/l"
    elif sbs_var in {"oxsPS"}:
        kvar_name = "oxsPS"
        kvar_units = r"% saturation"
    elif sbs_var in {"oxsMm/Kg"}:
        kvar_name = "oxsMm_Kg"
        kvar_units = "umol/kg"
    elif sbs_var in {"oxsdOC/dT"}:
        kvar_name = "oxsdOC_dT"
        kvar_units = "doc/dt"
    # IOW OXYGEN SATURATIONS
    elif sbs_var in {"iowOxML/L"}:
        kvar_name = "iowOxML/L"
        kvar_units = "ml/l"
    # GARCIA & GORDON OXYGEN SATURATION
    elif sbs_var in {"oxsolML/L"}:
        kvar_name = "oxsolML_L"
        kvar_units = "ml/l"
    elif sbs_var in {"oxsolMg/L"}:
        kvar_name = "oxsolMg_L"
        kvar_units = "mg/l"
    elif sbs_var in {"oxsolMm/Kg"}:
        kvar_name = "oxsolMm_Kg"
        kvar_units = "umol/kg"
    # WEISS OXYGEN SATURATION
    elif sbs_var in {"oxsatML/L"}:
        kvar_name = "oxsatML_L"
        kvar_units = "ml/l"
    elif sbs_var in {"oxsatMg/L"}:
        kvar_name = "oxsatMg_L"
        kvar_units = "mg/l"
    elif sbs_var in {"oxsatMm/Kg"}:
        kvar_name = "oxsatMm_Kg"
        kvar_units = "umol/kg"

    # OPTICAL SENSORS
    # TURNER CYCLOPS
    # CDOM
    elif sbs_var in {"cdomflTC0", "cdomflTC1", "cdomflTCdiff"}:
        kvar_name = sbs_var
        kvar_units = "ppb QS"
    # FLUORESCENCE
    elif sbs_var in {"chloroflTC0", "chloroflTC1", "chloroflTCdiff"}:
        kvar_name = sbs_var
        kvar_units = "ug/l"
    # CRUDE OIL
    elif sbs_var in {"croilflTC0", "croilflTC1", "croilflTCdiff"}:
        kvar_name = sbs_var
        kvar_units = "ppb QS"
    # FLUORESCEIN
    elif sbs_var in {"flflTC0", "flflTC1", "flflTCdiff"}:
        kvar_name = sbs_var
        kvar_units = "ppb"
    # OPTICAL BRIGHTENERS
    elif sbs_var in {"obrflTC0", "obrflTC1", "obrflTCdiff"}:
        kvar_name = sbs_var
        kvar_units = "ppb QS"
    # PHYCOCYANIN
    elif sbs_var in {"phycyflTC0", "phycyflTC1", "phycyflTCdiff"}:
        kvar_name = sbs_var
        kvar_units = "RFU"
    # PHYCOERYTHRIN
    elif sbs_var in {"phyeryflTC0", "phyeryflTC1", "phyeryflTCdiff"}:
        kvar_name = sbs_var
        kvar_units = "RFU"
    # REFINED FUELS
    elif sbs_var in {"rfuels0", "rfuels1", "rfuelsdiff"}:
        kvar_name = sbs_var
        kvar_units = "ppb NS"
    # RHODAMINE
    elif sbs_var in {"rhodflTC0", "rhodflTC1", "rhodflTCdiff"}:
        kvar_name = sbs_var
        kvar_units = "ppb"
    # TURNER CYCLOPS TURBIDITY
    elif sbs_var in {"turbflTC0", "turbflTC1", "turbflTCdiff"}:
        kvar_name = sbs_var
        kvar_units = "NTU"

    # BIOSPHERICAL
    # FLUORESCENCE
    elif sbs_var in {"chConctr", "naFluor", "product"}:
        kvar_name = sbs_var
        kvar_units = ""
    # CHELSEA AQUA 3
    # FLUORESCENCE
    elif sbs_var in {"flC", "flC1", "flCdiff"}:
        kvar_name = sbs_var
        kvar_units = "ug/l"
    # CHELSEA MINI CHL CON
    # FLUORESCENCE
    elif sbs_var in {"flCM"}:
        kvar_name = "flCM"
        kvar_units = "ug/l"
    # CHELSEA UV AQUATRACKA
    # FLUORESCENCE
    elif sbs_var in {"flCUVA", "flCUVA1", "flCUVAdiff"}:
        kvar_name = sbs_var
        kvar_units = "ug/l"
    # DR HAARDT
    # Fluorescence: Chlorophyll a, Phycoerythrin, Yellow Sub
    elif sbs_var in {"haardtC", "haardtP", "haardtY"}:
        kvar_name = sbs_var
        kvar_units = ""
    # SEAPOINT
    # FLUORESCENCE
    elif sbs_var in {"flSP", "flSP1", "flSPdiff"}:
        kvar_name = sbs_var
        kvar_units = ""
    # RHODAMINE
    elif sbs_var in {"flSPR"}:
        kvar_name = sbs_var
        kvar_units = ""
    # ULTRAVIOLET
    elif sbs_var in {"flSPuv0", "flSPuv1", "flSPuvdiff"}:
        kvar_name = sbs_var
        kvar_units = ""
    # SEATECH FLUORESCENCE
    elif sbs_var in {"flS"}:
        kvar_name = sbs_var
        kvar_units = ""
    # TURNER 10-005
    elif sbs_var in {"flT"}:
        kvar_name = "flT"
        kvar_units = ""
    # TURNER 10-Au-005
    elif sbs_var in {"flTAu"}:
        kvar_name = "flTAu"
        kvar_units = ""
    # TURNER SCUFA CORRECTED
    elif sbs_var in {"flSCC", "flSCC1", "flSCCdiff"}:
        kvar_name = sbs_var
        kvar_units = "RFU"
    # TURNER SCUFA
    elif sbs_var in {"flScufa", "flScufa1", "flScufadiff"}:
        kvar_name = sbs_var
        kvar_units = ""
    # WETLABS AC3
    elif sbs_var in {"wetChAbs"}:
        kvar_name = "wetChAbs"
        kvar_units = "1/m"
    # WETLABS CDOM
    elif sbs_var in {
        "wetCDOM",
        "wetCDOM1",
        "wetCDOM2",
        "wetCDOM3",
        "wetCDOM4",
        "wetCDOM5",
        "wetCDOMdiff",
    }:
        kvar_name = sbs_var
        kvar_units = "mg/m^3"
    # WETLABS AC3 CHLOROPHYLL
    elif sbs_var in {"wetChConc"}:
        kvar_name = "wetChConc"
        kvar_units = "mg/m^3"
    # WETLABS ECO-AFL
    elif sbs_var in {"flECO-AFL"}:
        kvar_name = "flECO_AFL"
        kvar_units = "mg/m^3"
    elif sbs_var in {"flECO-AFL1"}:
        kvar_name = "flECO_AFL1"
        kvar_units = "mg/m^3"
    elif sbs_var in {"flECO-AFL2"}:
        kvar_name = "flECO_AFL2"
        kvar_units = "mg/m^3"
    elif sbs_var in {"flECO-AFL3"}:
        kvar_name = "flECO_AFL3"
        kvar_units = "mg/m^3"
    elif sbs_var in {"flECO-AFL4"}:
        kvar_name = "flECO_AFL4"
        kvar_units = "mg/m^3"
    elif sbs_var in {"flECO-AFL5"}:
        kvar_name = "flECO_AFL5"
        kvar_units = "mg/m^3"
    elif sbs_var in {"flECO-AFLdiff"}:
        kvar_name = "flECO_AFLdiff"
        kvar_units = "mg/m^3"
    # WETLABS SEAOWL
    # FLUOROMETER
    elif sbs_var in {"flWETSeaOWLchl0", "flWETSeaOWLchl1", "flWETSeaOWLchldiff"}:
        kvar_name = sbs_var
        kvar_units = "?g/l"
    # FDOM
    elif sbs_var in {"flWETSeaOWLfdom0", "flWETSeaOWLfdom1", "flWETSeaOWLfdomdiff"}:
        kvar_name = sbs_var
        kvar_units = "?g/l"
    # WETSTAR
    elif sbs_var in {
        "wetStar",
        "wetStar1",
        "wetStar2",
        "wetStar3",
        "wetStar4",
        "wetStar5",
        "wetStardiff",
    }:
        kvar_name = sbs_var
        kvar_units = "mg/m^3"

    # BACKSCATTER
    # D & A BACKSCATTERANCE
    elif sbs_var in {"obs", "obs1", "obsdiff"}:
        kvar_name = sbs_var
        kvar_units = "NTU"
    # CHELSEA NEPHELOMETER
    elif sbs_var in {"nephc"}:
        kvar_name = "nephc"
        kvar_units = "FTU"
    # D & A
    elif sbs_var in {"obs3+"}:
        kvar_name = "obs3plus"
        kvar_units = "NTU"
    elif sbs_var in {"obs3+1"}:
        kvar_name = "obs3plus1"
        kvar_units = "NTU"
    elif sbs_var in {"obs3+diff"}:
        kvar_name = "obs3plusdiff"
        kvar_units = "NTU"
    # DR. HAARDT
    elif sbs_var in {"haardtT"}:
        kvar_name = "haardtT"
        kvar_units = ""
    # IFREMER
    elif sbs_var in {"diff"}:
        kvar_name = "diff"
        kvar_units = ""
    # SEATECH LS6000
    elif sbs_var in {"stLs6000", "stLs60001", "stLs6000diff"}:
        kvar_name = sbs_var
        kvar_units = ""
    # TURNER SCUFA
    elif sbs_var in {"obsscufa", "obsscufa1", "obsscufadiff"}:
        kvar_name = sbs_var
        kvar_units = ""
    # WETLABS RAW COUNTS
    elif sbs_var in {"wl0", "wl1", "wl2", "wl3", "wl4", "wl5"}:
        kvar_name = sbs_var
        kvar_units = "Counts"

    # TURBIDITY
    # SEAPOINT TURBIDITY
    elif sbs_var in {"seaTurbMtr", "seaTurbMtr1", "seaTurbMtrdiff"}:
        kvar_name = sbs_var
        kvar_units = "FTU"
    # WETLABS ECO BB
    elif sbs_var in {
        "turbWETbb0",
        "turbWETbb1",
        "turbWETbb2",
        "turbWETbb3",
        "turbWETbb4",
        "turbWETbb5",
        "turbWETbbdiff",
    }:
        kvar_name = sbs_var
        kvar_units = "m^-1/sr"
    # WETLABS ECO
    elif sbs_var in {
        "turbWETntu0",
        "turbWETntu1",
        "turbWETntu2",
        "turbWETntu3",
        "turbWETntu4",
        "turbWETntu5",
        "turbWETntudiff",
    }:
        kvar_name = sbs_var
        kvar_units = "NTU"
    # WETLABS SEAOWL
    elif sbs_var in {"turbWETSeaOWLbb0", "turbWETSeaOWLbb1", "turbWETSeaOWLbbdiff"}:
        kvar_name = sbs_var
        kvar_units = "m^-1/sr"

    # TRANSMISSOMETERS
    # CHELSEA/SEATECH BEAM ATTENUATION
    elif sbs_var in {"bat", "bat1", "batdiff"}:
        kvar_name = sbs_var
        kvar_units = "1/m"
    # WETLABS AC3 BEAM TRANSMISSION
    elif sbs_var in {"wetBAttn"}:
        kvar_name = "wetBAttn"
        kvar_units = "1/m"
    # WETLABS C-STAR BEAM ATTENUATION
    elif sbs_var in {
        "CStarAt",
        "CStarAt0",
        "CStarAt1",
        "CStarAt2",
        "CStarAt3",
        "CStarAt4",
        "CStarAt5",
        "CStarAtdiff",
    }:
        kvar_name = sbs_var
        kvar_units = "1/m"
    # CHELSEA/SEATECH BEAM TRANSMISSION
    elif sbs_var in {"xmiss", "xmiss1", "xmissdiff"}:
        kvar_name = sbs_var
        kvar_units = r"%"
    # WETLABS AC3 BEAM TRANSMISSION
    elif sbs_var in {"wetBTrans"}:
        kvar_name = "wetBTrans"
        kvar_units = r"%"
    # WETLABS C-STAR BEAM TRANSMISSION
    elif sbs_var in {
        "CStarTr",
        "CStarTr0",
        "CStarTr1",
        "CStarTr2",
        "CStarTr3",
        "CStarTr4",
        "CStarTr5",
        "CStarTrdiff",
    }:
        kvar_name = sbs_var
        kvar_units = r"%"

    # LISST
    elif sbs_var in {"lisstBC"}:
        kvar_name = "lisstBC"
        kvar_units = "1/m"
    elif sbs_var in {"lisstOT"}:
        kvar_name = "lisstOT"
        kvar_units = r"%"
    elif sbs_var in {"lisstMD"}:
        kvar_name = "lisstMD"
        kvar_units = "u"
    elif sbs_var in {"lisstTVC"}:
        kvar_name = "lisstTVC"
        kvar_units = "ul/l"

    # PAR
    elif sbs_var in {"cpar"}:
        kvar_name = "cpar"
        kvar_units = r"%"
    elif sbs_var in {"par"}:
        kvar_name = "par"
        kvar_units = ""
    elif sbs_var in {"par1"}:
        kvar_name = "par1"
        kvar_units = ""
    elif sbs_var in {"par/log"}:
        kvar_name = "par_log"
        kvar_units = "umol photons/m2/s"
    elif sbs_var in {"spar"}:
        kvar_name = "spar"
        kvar_units = ""

    # BIOGEOCHEMICAL SENSORS
    # METHANE FRANATECH METS
    elif sbs_var in {"meth"}:
        kvar_name = "meth"
        kvar_units = "umol/l"
    elif sbs_var in {"methT"}:
        kvar_name = "methT"
        kvar_units = "deg C"
    # NITROGEN
    elif sbs_var in {"n2satML/L"}:
        kvar_name = "n2satML_L"
        kvar_units = "ml/l"
    elif sbs_var in {"n2satMg/L"}:
        kvar_name = "n2satMg_L"
        kvar_units = "mg/l"
    elif sbs_var in {"n2satumol/kg"}:
        kvar_name = "n2satumol_kg"
        kvar_units = "umol/kg"
    # OXIDATION REDUCTION POTENTIAL
    elif sbs_var in {"orp"}:
        kvar_name = "orp"
        kvar_units = "mV"
    # pH
    elif sbs_var in {"ph", "phInt", "phExt"}:
        kvar_name = sbs_var
        kvar_units = ""
    # ZAPS
    elif sbs_var in {"zaps"}:
        kvar_name = "zaps"
        kvar_units = "nmol"

    # USER DEFINED VARIABLE
    elif sbs_var in {"user", "user1", "user2", "user3", "user4", "user5"}:
        kvar_name = sbs_var
        kvar_units = ""
    # USER POLYNOMIAL
    elif sbs_var in {"upoly0", "upoly1", "upoly2"}:
        kvar_name = sbs_var
        kvar_units = ""

    # PACKAGE STATUS
    # ALTIMETER
    elif sbs_var in {"altM"}:
        kvar_name = "altM"
        kvar_units = "m"
    elif sbs_var in {"altF"}:
        kvar_name = "altF"
        kvar_units = "ft"

        # Bottom Contact
    elif sbs_var in {"bct"}:
        kvar_name = "bct"
        kvar_units = ""

    # WATER SAMPLER STATUS
    # BOTTLE POSITION, BOTTLES CLOSED (HB), BOTTLES FIRED, NEW POSITION'
    elif sbs_var in {"bpos", "HBBotCls", "nbf", "newpos"}:
        kvar_name = sbs_var
        kvar_units = ""
    # PUMP STATUS
    elif sbs_var in {"pumps"}:
        kvar_name = "pumps"
        kvar_units = ""

    # DATA STREAM
    # BYTE COUNT
    elif sbs_var in {"nbytes"}:
        kvar_name = "nbytes"
        kvar_units = ""
    # MODULO ERROR COUNT
    elif sbs_var in {"modError"}:
        kvar_name = "modError"
        kvar_units = ""
    # MODULO WORD
    elif sbs_var in {"mod"}:
        kvar_name = "mod"
        kvar_format = r"%s"
        kvar_units = ""
    # SCAN COUNT
    elif sbs_var in {"scan"}:
        kvar_name = "scan"
        kvar_units = ""
    # SCANS PER BIN
    elif sbs_var in {"nbin"}:
        kvar_name = "nbin"
        kvar_units = ""

    # FLAGS
    elif sbs_var in {"flag"}:
        kvar_name = "flag"
        kvar_units = ""

    else:
        # assign the original SBS name to the output variable
        kvar_name = sbs_var
        # find and replace all the characters that are not alphabetic,
        # numbers or underscores with underscores
        kvar_name = re.sub(r"\W", "_", kvar_name)
        print(f"Assigning {sbs_var} as {kvar_name}")
        kvar_units = ""

    #     elif sbs_var in  ''
    #         kvar_name = ''
    # #         kvar_units = ''

    return {"name": kvar_name, "format": kvar_format, "units": kvar_units}
