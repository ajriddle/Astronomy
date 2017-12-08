# Work-Code
This repo is where you can find code and results from my time as an astronomy grad student at UT.

## Results Files
There are several pairs of .png files that show some example results of my work on eclipsing binary stars. Images with 'triangle' in the name are the results of a Markov Chain Monte Carlo analytical fit to observed data (triangle plots).

<img src=./Results/MG1-2062466_triangle.png alt="triangle plot" width=400 height=400> 

The companion files to the triangle plots (same name without the "triangle" part) show the fitted model plotted over the observed radial velocity data. The bottom panel of this figure shows the residuals of the model along with the observed errors of the data (the errorbars). The different colors represent the two stars of the binary system.
<img src=./Results/MG1-2062466.png alt="RV plot" width=500 height=400>

## Code
- [RV_MCMC.py](./Code/RV_MCMC) - Performs MCMC analysis to fit 5-7 parameter orbital model to observed radial velocity data.
- [eclipse.py](./Code/eclipse.py) - analyzes light curve data taken from the 0.8-m telescope at McDonald Observatory, looking for stellar eclipses
- [yygem_phases.py](./Code/yygem_phases.py) - creates stellar phases for the binary system YY Gem from historical data
- [todcor.py](./Code/todcor.py) - main function to compute the two-dimensional cross-correlation (to determine radial velocities) from an observed spectrum of a binary star. It employs the following functions that handle data from different sources to do the bulk of the actual calculations:
    - [cps_fcns.py](./Code/cps_fcns.py)
    - [archive_fcns.py](./Code/archive_fcns.py)
- [shift_estimate.py](./Code/shift_estimate.py) - Loops through spectroscopic observations of a system one epoch at a time and then through different spectral orders for each epoch, showing the 1-D cross-correlation for each order to give an estimate of the observed radial velocities.
- [order_scatter.py](./Code/order_scatter.py) - Calculates radial velocity scatter for different spectral orders. This was used to determine which orders produced the best radial velocity measurements.
- [masses.py](./Code/masses.py) - Calculates stellar component masses for a binary system based on the parameters found from the MCMC analysis.
- [fits_read.py](./Code/fits_read.py) - Reads in .fits image files and converts them to numpy arrays that contain flux versus wavelength data for the different orders (rows).
- [EB_veloity.py](./Code/EB_velocity.py) Converts from stellar orbital parameters to radial velocities, which invloves solving a transcendental equation numericaly.
- [read_rvs.py](./Code/read_rvs.py) - Reads in radial velocities for one system and performs linear regression to fit a simple orbital model. This only worked for systems with known eclipse timing as linear regression would only work on those systems.
 - [phases.py](./Code/phases.py) - Predicts eclipse phases for a set of binary systems for a defined observing run.
 - [phase_calculate.py](./Code/phase_calculate.py) - Extracts/converts time of minimum (eclipse time) and orbital period from different sources to a standardized format. 