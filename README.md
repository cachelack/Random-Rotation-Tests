# Random-Rotation-Tests
A collection of R code for performing statistical tests using random rotations from SO(n)

## quadFormTest.r
This is a collection of R functions for performing a random rotation test for testing
for autocorrelation in univariate time series data.  In particular the function
doSOTest() returns a p-value for testing a time series (dat) at lag (lag) for autocorrelation.
The other functions included perform a variety of simulations and compare our doSOTest() function
against the classic Durbin-Watson and Breusch-Godfrey tests.

## observed-solar-cycle-indices.json
A time series of solar sunspots and intensity publicly available at https://www.swpc.noaa.gov/products/solar-cycle-progression
but included here for convenience.
