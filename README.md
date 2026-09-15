# MATLAB Code for Perturbed Weddle Quadrature Numerical Experiments

This repository contains the MATLAB implementation associated with the manuscript:

**A New Error Estimation Framework for Perturbed Weddle Quadrature Using Second Derivatives**

The code is provided to improve the reproducibility of the numerical experiments presented in the manuscript.

## Numerical Methods

The MATLAB implementation includes the following quadrature methods:

- NC-2: Composite three-node Newton-Cotes (Simpson's) rule
- G-2: Composite two-point Gauss-Legendre quadrature
- M-2: Mixed Newton-Cotes-Gaussian quadrature method

## Test Problems

The numerical experiments include:

1. Polynomial function
2. Exponential function
3. Trigonometric function
4. Counterexample
5. q-Digamma function
6. Modified Bessel function

## Numerical Resolution

The computations are performed using:

N = 12, 24, 48, 96, 192, 384.

The code evaluates the absolute quadrature errors and observed convergence
orders for the numerical methods considered in the manuscript.

## Reproducibility

Running `APNUM_Appendix_NumericalCode.m` reproduces the numerical data and
convergence results reported in the revised manuscript.

## Software

The computations are implemented in MATLAB.
