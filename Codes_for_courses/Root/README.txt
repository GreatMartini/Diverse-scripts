This code is a "ROOT" (c++ based data analysis framework) code that finds the best estimator of a parameter for a signal embedded in noise. This code is part of a project made for a course on signal detection.

First, the "construct()" function creates a gaussian signal embedded in a decreasing exponential noise function and, transforms the generated signal into a histogram.

The "likelihood_imp(double k)" implements the poissonian likelihood function of the signal for the parameter k in the model histogram to adjust it to the signal.

The "minimisation()" function retrieves the likelihood function and finds the minimum (in this case maximum because of the implmentation) of the parabolic shape of the likelihood function in order to minimize it and find the best estimator for the parameter k. Afterwards it calculates the pulls and draws a histogram.