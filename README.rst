SED Fitting Tutorial
====================

This repository contains a series of iPython notebooks which will teach you how to perform very basic SED fitting. The assumption is that you're familiar with the basic concepts (e.g. by having read Section 2 of `Conroy 2013 <https://arxiv.org/abs/1301.7095>`_), and now want to start implementing those concepts in Python.

I wrote it in the hope that it might be useful either for students who wanted to start writing their own SED fitting package, or for people who are using other people's SED fitting tools and want to understand a little better how they work.

The method the code applies is grid-based chi-squared minimisation. It fits age, mass and redshift for simple stellar population models from BPASS. The code is deliberately oversimplified and inefficient in order to make the syntax as clear as possible. 

If you're working off this code as a starting point, suggestions for developing it further are as follows:
 - Vectorise mathematical operations with numpy to speed up the code, although beware running out of memory.
 - Perform analytic chi-squared minimisation to determine the best fitting mass and only grid over age and redshift.
 - Add a function to apply the Calzetti et al. (2000) dust attenuation law to the spectra and add Av as another parameter.
 - Add the ability to fit different metallicities and more complex star-formation histories by adding SSP models together.
 - Put the whole thing inside a class so you're not accessing variables from the global namespace within functions.
 - Replace the grid search method with some sort of functional minimisation (scipy.optimise) or MCMC (emcee) routine.

At that point you've already got a highly competitive SED fitting code, where you go from there is up to you...

Or you could just used `Bagpipes <https://github.com/ACCarnall/bagpipes>`_.