//
//  ChockSample.hpp
//  AvoidReject
//
//  Created by Mitchell, Scott A on 7/28/21.
//  Copyright © 2021 Mitchell, Scott A. All rights reserved.
//

#ifndef ChockSample_hpp
#define ChockSample_hpp

#include <cmath>
#include <assert.h>
#include <limits>

// random number generation doesn't belong here, assume u random variable is passed in

// should be defined in cmath or math.h
#ifndef M_PI
#define M_PI        3.14159265358979323846264338327950288   /* pi             */
#endif
#ifndef M_PI_2
#define M_PI_2      1.57079632679489661923132169163975144   /* pi/2           */
#endif
#ifndef M_PI_4
#define M_PI_4      0.785398163397448309615660845819875721  /* pi/4           */
#endif
#ifndef M_1_PI
#define M_1_PI      0.318309886183790671537767526745028724  /* 1/pi           */
#endif

using std::tan;
using std::cos;

// ========== interface

// sample phi and radius uniformly by area, given uniform random variables u_phi and u_r
void sample_chock( double u_phi, double u_r, double phi, double &phi_s, double &r_s);

// find phi_sample : area_sample = u*area_total(phi)
// Inverse CDF transform sampling
void sample_phi(double u, double phi_total, double &phi_sample, double &area_total, double &area_sample);

// find r_sample : swept_area(r_s) = u*swept_area( r(phi))
// Inverse CDF transform sampling
void sample_r( double u, double phi, double &r_s);


// ========== implementation

inline
void sample_r( double u, double phi, double &r_s)
{
  assert( fabs(phi) + 2.*std::numeric_limits<double>::epsilon() < M_PI_2);
  assert( u >= 0. );
  assert( u <= 1. );
  // double &r_total = chock_r(phi) // not needed
  const double tanphi = tan(phi);
  r_s = sqrt( u * tanphi * tanphi + 1.);
}

inline
void sample_chock( double u_phi, double u_r, double phi, double &phi_s, double &r_s)
{
  double A, As;
  sample_phi(u_phi, phi, phi_s, A, As);
  sample_r(u_r,phi_s,r_s);
}



#endif /* ChockSample_hpp */
