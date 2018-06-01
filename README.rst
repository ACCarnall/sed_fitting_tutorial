SED Fitting Tutorial
====================

This repository contains a series of iPython notebooks which will teach you how to perform very basic SED fitting.

I wrote it in the hope that it might be useful for students who wanted to start writing their own SED fitting package.

The method the code applies is grid-based chi-squared minimisation. It fits age, mass and redshift for simple stellar population models from BPASS.

The code is deliberately oversimplified and inefficient in order to make the syntax as clear as possible. 

If you're working off this code as a starting point, suggestions for developing it further are as follows:
 - Vectorise mathematical operations with numpy to speed up the code, although beware running out of memory.
 - Perform analytic chi-squared minimisation to determine the best fitting mass and only grid over age and redshift. Again this should produce a big speed boost.
 - Add a function to apply the Calzetti et al. (2000) dust attenuation law to the spectra and add Av as another grid parameter, this will drastically improve the quality of your results.
 - Add the ability to fit different metallicity grids and more complex star-formation histories by adding SSP models together.
 - Put the whole thing inside a class so you're not accessing variables from the global namespace within functions.
 - Replace the grid search method with some sort of functional minimisation (scipy.optimise) or MCMC (emcee) routine, again for speed.

At that point you've already got a highly competitive SED fitting code, where you go from there is up to you...