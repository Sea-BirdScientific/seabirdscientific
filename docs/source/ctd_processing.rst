.. _ctd_processing:

CTD Processing
##############

Low Pass Filter
***************

The low-pass filter removes high frequency data to make the data smoother. To make sure there is zero phase shift (no time shift), the filter is first applied forward through the data and then backward through the forward-filtered data. This removes any time shifts caused by the filter.  

Pressure data is typically filtered with a time constant equal to four times the CTD scan rate. Conductivity and temperature are typically filtered for *some* CTDs. Typical time constants are:  

=======================  ===========  ============  ========
Instrument               Temperature  Conductivity  Pressure
=======================  ===========  ============  ========
SBE 19plus or 19plus V2  0.5s         0.5s          1.0s    
=======================  ===========  ============  ========

Filter Formula
==============

For a low-pass filter with time constant :math:`Γ`:  

    | :math:`Γ= 1/2πf`

    | :math:`T = sample\ interval\ (seconds)`

    | :math:`S_0 = 1/Γ`

Laplace transform of the transfer function of a low-pass filter (single pole) with a time constant of :math:`Γ` seconds is:  

    | :math:`H(s) = \frac{1}{1 + (S/S_0)}`
  
Using the bilinear transform:  

    | :math:`S - f(z) \overset{\Delta}{=} \frac{2(1-z^{-1})}{T(1+z^{-1})} = \frac{2(z-1)}{T(z+1)}`

    | :math:`H(z) = \frac{1}{1 + \frac{2(z-1)}{T(z+1)S_0}} = \frac{z^{-1}+1}{1 + \frac{2}{TS_0}(1+\frac{1-2/TS_0}{1+2/TS_0}z^{-1})}`

If:  

    | :math:`A = \frac{1}{1 + \frac{2}{TS_0}}`

    | :math:`B = \frac{1 - \frac{2}{TS_0}}{1 + \frac{2}{TS_0}}`

Then:  

    | :math:`H(z) = \frac{Y(z)}{X(z)} = \frac{A(z^{-1}+1)}{1+Bz^{-1}}`
  
Where :math:`z^{-1}` is the unit delay (one scan behind).  

    | :math:`Y(z) (1 + Bz -1 ) = X(z) A (z -1 + 1)`

    | :math:`y[N] + By[N-1] = Ax[N-1] + Ax[N]`

    | :math:`y[N] = A(x[N] + x[N-1]) - By[N-1]`

****

Align CTD
*********

The Align CTD function aligns parameter data in time, relative to pressure, so that the calculated salinity, dissolved oxygen concentration, and other parameters are made with measurements from the same parcel of water. Typically, Align CTD is used to align temperature, conductivity, and oxygen measurements relative to pressure.  

There are three principal causes of misalignment of CTD measurements:

- physical misalignment of the sensors in depth
- inherent time delay (time constants) of the sensor responses
- water transit time delay in the pumped plumbing line - the time it takes the parcel of water to go through the plumbing to each sensor (or, for free-flushing sensors, the related flushing delay, which depends on profiling speed)  

When measurements are correctly aligned, salinity spikes (and density) errors are minimized, and oxygen data agrees with the correct pressure (e.g., temperature vs. oxygen plots agree between down and up profiles).  

Conductivity and Temperature
============================

Temperature and conductivity are often misaligned with respect to pressure. It can help to move the temperature and conductivity relative to pressure.  

The best diagnostic of correct alignment is the removal of salinity spikes that align with very sharp temperature steps. To determine the best alignment, make a plot of 10 meters of temperature and salinity data at a depth that contains a very sharp temperature step. For the downcast, when temperature and salinity decrease with increased pressure:

- A negative salinity spike at the conductivity step means that conductivity leads temperature (conductivity sensor sees step before temperature sensor does). Move conductivity relative to temperature a negative number of seconds.  
- If the salinity spike is positive, move conductivity relative to temperature a positive number of seconds.

The best alignment of conductivity and temperature occurs when the salinity spikes are minimized. It may be necessary to try different intervals to find the best alignment.

Typical Temperature Alignment
-----------------------------

The SBE 19, 19plus, and 19plus V2 use a temperature sensor with a 0.5 second response time, which is slower than the conductivity sensor.

===================  ===========================================
Instrument           Advance of temperature relative to pressure
===================  ===========================================
19plus or 19plus V2  +0.5 seconds
===================  ===========================================


Typical Conductivity Alignment
------------------------------

For an SBE 19plus or 19plus V2 with a standard 2000-rpm pump, do not move conductivity.  
However, if temperature is moved relative to pressure and you do not want to change the relative timing of temperature and conductivity, you must add the same interval to conductivity.  

Oxygen
======

Oxygen data is also systematically delayed with respect to pressure. The two primary causes are the long time constant of the oxygen sensor (for the SBE 43, ranging from 2 seconds at 25 ºC to approximately 5 seconds at 0 ºC) and an additional delay from the transit time of water in the pumped plumbing line. As with temperature and conductivity, you can move oxygen data relative to pressure to adjust for this delay.  

===================  ======================================
Instrument           Advance of oxygen relative to pressure
===================  ======================================
19plus or 19plus V2  +3 to +7 seconds
===================  ======================================

****

Cell Thermal Mass
*****************

The Cell Thermal Mass function uses a recursive filter to remove conductivity cell thermal mass effects from the measured conductivity. Typical values for alpha and 1/beta are:  

======================================================  =====  ======
Instrument                                              alpha  1/beta
======================================================  =====  ======
SBE 19plus or 19plus V2 with TC duct and 2000 rpm pump  0.04   8.0
======================================================  =====  ======

Cell Thermal Mass Formulas
==========================

    | :math:`a = \frac{2*alpha}{(sample\ interval * beta + 2)}`

    | :math:`b = 1 - \frac{2a}{alpha}`

    | :math:`\frac{dc}{dT} = 0.1 + 6e\text{-}4*(temperature-20)`

    | :math:`dT = temperature - previous\ temperature`

    | :math:`ctm = -b*previous\ ctm+a*\frac{dc}{dT}*dT`

Where sample interval is in seconds, temperature is in °C, ctm is in S/m, amplitude is the alpha value, and time_constant is 1/beta.

To determine the values for alpha and beta, see:

- The Correction for Thermal-Lag Effects in Sea-Bird Scientific CTD Data (https://journals.ametsoc.org/view/journals/atot/11/4/1520-0426_1994_011_1151_tcftle_2_0_co_2.xml), Morison, J., R. Andersen, Larson, N., D'Asaro, E., and Boyd, T.,  Journal of Atmospheric and Oceanic Technology (JAOT), V11(4), August 1994, 1151-1164
- Thermal Inertia of Conductivity Cells: Theory (https://journals.ametsoc.org/view/journals/atot/7/5/1520-0426_1990_007_0741_tiocct_2_0_co_2.xml), Lueck, R.G., Journal of Atmospheric and Oceanic Technology (JAOT), V7(5), Oct 1990, 741-755.

****

Loop Edit
*********

The Loop Edit function identifies scans as bad. The flag value associated with the scan is set to badflag in input .cnv files that have pressure slowdowns or reversals (typically caused by ship heave). Optionally, Loop Edit can also identify scans related to an initial surface soak with badflag. The badflag value defaults to -9.99e-29.

Loop Edit operates on three successive scans to determine velocity. This is such a fine scale that noise in the pressure channel from counting jitter or other unknown sources can cause Loop Edit to mark scans with badflag in error. Therefore, you must do the Filter on the pressure data to reduce noise before you do the Loop Edit. See Filter for pressure filter recommendations for each instrument.

****

Derive TEOS-10
**************

The Derive TEOS-10 function uses temperature, conductivity, practical salinity, pressure, latitude, and/or longitude to compute the following thermodynamic parameters using TEOS-10 equations:  

The table below references python functions from `TEOS-10 <https://github.com/TEOS-10/GSW-Python>`_ that are the same functions in SBE Data Processing. Note: Results may differ slightly from SBE Data Processing due to GSW version differences.

==================================  ================================
Variable Name                       GSW-Python function             
==================================  ================================
Absolute Salinity                   gsw.SA_from_SP                  
Absolute Salinity Anomaly           gsw.deltaSA_from_SP             
Adiabatic Lapse Rate                gsw.adiabatic_lapse_rate_from_CT
Conservative Temperature            gsw.CT_from_t                   
Conservative Temperature, Freezing  gsw.CT_freezing                 
Density, TEOS-10                    gsw.rho \*                      
Dynamic Enthalpy                    gsw.dynamic_enthalpy            
Enthalpy                            gsw.enthalpy                    
Entropy                             gsw.entropy_from_t              
Gravity                             gsw.grav                        
Internal Energy                     gsw.internal_energy             
Isentropic Compressibility          gsw.kappa                       
Latent Heat of Evaporation          gsw.latentheat_evap_CT          
Latent Heat of Melting              gsw.latentheat_melting          
Potential Temperature               gsw.pt0_from_t                  
Preformed Salinity                  gsw.Sstar_from_SA               
Reference Salinity                  gsw.SR_from_SP                  
Saline Contraction Coefficient      gsw.beta                        
Sound Speed                         gsw.sound_speed                 
Specific Volume                     gsw.specvol                     
Specific Volume Anomaly             gsw.specvol_anom_standard       
Temperature Freezing                gsw.t_freezing                  
Thermal Expansion Coefficient       gsw.alpha                       
==================================  ================================

\* use gsw.rho with reference pressure for the sigmas

****

Wild Edit
*********

The Wild Edit function identifies wild points in the data. The data value is replaced with a badflag value, by default -9.99e-29. Wild Edit's algorithm requires two passes through the data. The first pass gives an accurate estimate of the data's true standard deviation, while the second pass replaces the appropriate data with badflag.

Wild Edit operates as follows:

1. Compute the mean and standard deviation of data in block (specified by Scans per Block) for each selected variable. Temporarily flag values that differ from mean by more than standard deviations specified for pass 1.
2. Recompute mean and standard deviation, and do not include temporarily flagged values. Identify the values that differ from the mean by more than standard deviations specified for pass 2 and replace the data value with badflag.
3. Repeat Steps 1 and 2 for next block of scans. If the last block of data in the input file has less than specified number of scans, use data from previous block to fill in the block.

If the data file is particularly corrupted, you may need to do Wild Edit more than once, with different block sizes and number of standard deviations. If the input file has some variables with large values and some with relatively smaller values, it may be necessary to run Wild Edit more than once, and change the value for Keep data within this distance of mean so that it is meaningful for each variable. Increase Scans per block from 100 to approximately 500 to get better results.

****

Window Filter
*************

The Window Filter function gives four types of window filters and a median filter that can be used to make the data smoother:  

- Window filters calculate a weighted average of data values about a center point and replace the data value at the center point with this average.  
- The median filter calculates a median for data values around a center point and replaces the data value at the center point with the median 

Descriptions and Formulas
=========================

Shape and length define filter windows:  

- Window Filter gives four window shapes: boxcar, cosine, triangle, and Gaussian.  
- The minimum window length is 1 scan, and the maximum is 511 scans. Window length must be an odd number, so that the window has a center point. If a window length is specified as an even number, Window Filter automatically adds 1 to make the length odd.  

The window filter calculates a weighted average of data values around a center point, using the transfer function below:  

    | :math:`y(n) = \displaystyle\sum_{k=-L/2}^{L/2}w(k)x(n-k)`

The window filter process is similar for all filter types:

1. Filter weights are calculated (see the equations below).
2. Filter weights are normalized to sum to 1.

   - When the data includes scans that have been identified as bad, and the ``exclude_flags`` variable is set to true, the weights are renormalized, and do not include the filter elements that would operate on the bad data point.  

In the equations below:  

    | :math:`L = window\ length\ in\ scans\ (must\ be\ odd)`

    | :math:`n = window\ index,\ -\frac{L-1}{2}..\frac{L-1}{2}`

    | :math:`w(n) = set\ of\ window\ weights`

Boxcar Filter
-------------
    | :math:`w(n) = \frac{1}{L}`

Cosine Filter
-------------
    | :math:`w(n) = cos(\frac{n\pi}{L+1})`

Triangle Filter
---------------
    | :math:`w(n) = 1-\frac{|2n|}{L+1}`

Gaussian Filter
---------------
    | :math:`phase = \frac{offset}{sample\ interval}`

    | :math:`scale = log(2)(2\frac{sample\ rate}{half\ width\ (scans)})^2`

    | :math:`w(n) = e^{(n-phase)^2 * scale}`

The Gaussian window has parameters of halfwidth (in scans) and offset (in time), in addition to window length (in scans). These extra parameters filter and move data in time in one operation. Halfwidth determines the width of the Gaussian curve. A window length of 9 and halfwidth of 4 gives a set of filter weights that fills the window. A window length of 17 and halfwidth of 4 gives a set of filter weights that fills only half the window. If the filter weights do not fill the window, the offset parameter may be used to move the weights within the window so that the edge of the Gaussian curve is not cut.  

Median Filter
-------------
The median filter does not make data smoother in the same sense as the window filters described above. The median filter is most useful in spike removal. A median value is determined for a specified window, and the data value at the window's center point is replaced by the median value.  