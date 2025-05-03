# ED_ST_GEOSTATS
Edinburgh Space Time Geostatistics

I wrote this programme many years ago as a tool for the completion of my PhD thesis. Ive since moved on to orher disciplines, but im sharingbit here in the hope it is useful to others.

The Fortran 90 routines allow for the spatio-temporal regionalization of data via geostatistical methods. The
programs utilize the product and product sum covariance representations of spatio-
temporal data interactions. The code described allows interpolation of a data set over
and arbitrarily spaced grid in continuous spatial and temporal coordinate systems. Error
analyses are provided via the Jack Knife. The resultant spatio-temporal fields represent
the expectation of a random function (RF), conditioned on the observed data and
covariance model. The techniques implemented allow production of fields of estimation
variance. We also provide code for simulation from the RF, via Sequential Gaussian
Simulation (SGS). The SGS technique makes random draws from the RE, and allows
the user to quantify the uncertainty of the interpolated field. SGS allows Monte Carlo
analysis of the regionalisation by ensuring draws from the distribution (described by the
expectation and estimation variance) conform to the observed spatio-temporal
covariance of the data. The SGS technique is particularly useful in cases where the user
intends to parameterise a model with a regionalised field, as the interpolation uncertainty
can be propagated through the model to produce appropriate confidence intervals.
Instructions for the use of the software are provided, along with sufficient background
theory to successfully implement spatio-temporal regionalisation