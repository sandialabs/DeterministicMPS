# DeterministicMPS
Maximal Poisson-disk Sampling in 2D, without approximations or rejections, in deterministic linear time, using chocks.  A chock is a curved region bounded by a disk, a tangent, and a ray.  Our formula for inverse transform sampling over chocks is numerically closed.

The software is deemed to be Publicly Available.  Copyright 2022 Sandia Corporation.  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.  SCR#:2749.0  

The main purpose of this software is to reproduce the results, data and figures, in the following publication:

"Deterministic Linear Time for Maximal Poisson-Disk Sampling using Chocks without Rejection or Approximation" by Scott A. Mitchell.  Conditionally accepted to Eurographics Symposium on Geometry Processing, SGP, 2022.


## Compiling
A CMakeLists.txt file is provided.  
The software has been compiled using Xcode on Mac OS X, and cmake and gcc on Linux.  The code uses c++11 extensions.  The code is self-contained and there are no library dependencies.  See the "main" method in the main.cpp file, which generates the figures (or data for the figures) from the paper.  There are other top-level methods, by default commented out, for unit testing and generating distributions over rectangles of different sizes and periodicity.
