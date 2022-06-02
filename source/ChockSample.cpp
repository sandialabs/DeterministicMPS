//
//  ChockSample.cpp
//  AvoidReject
//
//  Created by Mitchell, Scott A on 7/28/21.
//  Copyright © 2021 Mitchell, Scott A. All rights reserved.
//

#include "ChockSample.hpp"

#include <iostream>
#include <cmath>
#include <assert.h>

#include "Geometry.hpp"

using std::tan;
using std::pow;

void sample_phi(double u, double phi_total, double &phi_sample, double &area_total, double &area_sample)
{
  area_total = chock_area(phi_total);
  
  area_sample = u*area_total; // + std::numeric_limits<double>::epsilon();
  
  // explicitly check for zero to avoid divided-by-zero
  //   the method is numerically stable for area_sample = __DBL_EPSILON__
  if (area_sample == 0.)
  {
    phi_sample = 0.;
    // std::cout << "\ndebug: sampling zeroth-area of chock";
    // std::cout << "\nu " << u << ", phi_total " << phi_total << ", phi_sample " << phi_sample << ", area_total " << area_total << ", area_sample " << area_sample << std::endl;
    return;
  }
  
  // This is hard-coded for double precision to 5 iterations.
  //  
  // Higher precision requires more iterations,
  // but can be implemented to still cost only a constant times the cost of single std::tan() call.
  // In the first iteration use very low precision, say single precision.
  // In each iteration roughly double the precision, which doubles the cost as well.
  // Perform O(log b) iterations, where b is the bits of precision desired.
  // The sum of the cost of all these iterations is merely twice the cost of the last iteration,
  // through the sum of a geometric series.
  
  // phi 0 initial guess
  const auto p0=pow(6.*area_sample, 1./3.);
  
  // iter 1
  const auto t0=tan(p0);
  const auto f0=(t0-p0)/2. - area_sample;
  const auto fp0=t0*t0/2.; // derivative
  const auto p1=p0 - f0/fp0;
  
  // iter 2
  const auto t1=tan(p1);
  const auto f1=(t1-p1)/2. - area_sample;
  const auto fp1=t1*t1/2.;
  const auto p2=p1 - f1/fp1;
  
  // iter 3
  const auto t2=tan(p2);
  const auto f2=(t2-p2)/2. - area_sample;
  const auto fp2=t2*t2/2.;
  const auto p3=p2 - f2/fp2;
  
  // iter 4
  const auto t3=tan(p3);
  const auto f3=(t3-p3)/2. - area_sample;
  const auto fp3=t3*t3/2.;
  const auto p4=p3 - f3/fp3;
  
  // iter 5
  const auto t4=tan(p4);
  const auto f4=(t4-p4)/2. - area_sample;
  const auto fp4=t4*t4/2.;
  const auto p5=p4 - f4/fp4;
  
  // check that we converged
  // iter 6
  const auto t5=tan(p5);
  const auto f5=(t5-p5)/2. - area_sample;
  // const auto fp5=t5*t5/2.;
  // const auto p6=p5 - f5/fp5;
  
  phi_sample = p5;
  
  // debug
  if (/* DISABLES CODE */ (0))
  {
    std::cout << "\nu " << u << ", phi_total " << phi_total << ", phi_sample " << phi_sample << ", area_total " << area_total << ", area_sample " << area_sample << std::endl;
    std::cout << "p0 " << p0 << std::endl;
    std::cout << "t0 " << t0 << ", f0 " << f0 << ", fp0 " << fp0 << ", p1 " << p1 << std::endl;
    std::cout << "t1 " << t1 << ", f1 " << f1 << ", fp1 " << fp1 << ", p2 " << p2 << std::endl;
    std::cout << "t2 " << t2 << ", f2 " << f2 << ", fp2 " << fp2 << ", p3 " << p3 << std::endl;
    std::cout << "t3 " << t3 << ", f3 " << f3 << ", fp3 " << fp3 << ", p4 " << p4 << std::endl;
    std::cout << "t4 " << t4 << ", f4 " << f4 << ", fp4 " << fp4 << ", p5 " << p5 << std::endl;
    std::cout << "t5 " << t5 << ", f5 " << f5 << " should be 0 " << std::endl;
  }
  
  assert(fabs(f5)<=std::numeric_limits<double>::epsilon());
  assert(fabs(f5/phi_sample)<=2*std::numeric_limits<double>::epsilon()); // <2 is worst factor encountered so far
  
}
