So you got an idea for an EISCAT_3D experiment, and you want to estimate what sorts of uncertainties you'll get and what integration times are required?

You've come to the right place!

By default, the scripts in this repository will produce the synthetic EISCAT_3D data set described by Hatch et al (publication forthcoming). But they can easily be modified for different beam configurations, integration times, radar system parameters, etc.

The hitch is that you need to have both a working python setup _and_ a working R setup. (The R code runs a slightly modified version of I. Virtanen's ISgeometry package. This package is _necessary_ for estimating uncertainties of electron density, ion temperature, electron/ion temperature, and ion convection velocity.

## ´1_select_beam_geometry.py´

* Specify transmitter and receiver locations
* Pick number of beams, azimuths, elevations, and which altitudes to sample along each beam
* Output .txt files to be used by ´2_get_noise_estimates.R´

## ´2_get_noise_estimates.R´

* Read output from ´1_select_beam_geometry.py´
* Specify radar system parameters such as noise temperature, duty cycle, transmission frequency and power, and electron density and width for estimation of transmitter self-noise, etc.
* Output relative ACF noise levels

## ´3_get_gemini_for_all_timestamps_at_e3d_points.py´

* Coming soon!
* Big idea is to pull in synthetic/simulated data specified on some grid, and to 
** Resample these data at the points specified in ´1_select_beam_geometry.py´
** Estimate uncertainties by specifying an integration time and using ´ErrorTable_E3D_4params.txt´ as well as relative ACF noise levels produced by ´2_get_noise_estimates.R´
