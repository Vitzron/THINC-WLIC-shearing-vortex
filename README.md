# THINC-WLIC-shearing-vortex
shearing vortex simulation with THINC/WLIC

## Reference
Xiao F, Honma Y, Kono T. A simple algebraic interface capturing scheme using hyperbolic tangent function[J]. International Journal for Numerical Methods in Fluids, 2005, 48(9): 1023–1040.

Yokoi K. Efficient implementation of thinc scheme: a simple and practical smoothed vof algorithm[J]. Journal of Computational Physics, 2007, 226(2): 1985–2002.

## Eigen package REQUIRED!

See [Eigen](https://eigen.tuxfamily.org/) for more information about Eigen.

## make

Change the LIBS in makefile to the Eigen directory.

## ./main arg
arg: mesh number along each direction

## python plot.py [arg]
arg: figure name saved
