//
//  TestChockSample.cpp
//  AvoidReject
//
//  Created by Mitchell, Scott A on 7/28/21.
//  Copyright © 2021 Mitchell, Scott A. All rights reserved.
//

#include "TestChockSample.hpp"

#include <iostream>
#include <fstream>

#include "MyRand.hpp"
#include "ChockSample.hpp"
#include "Geometry.hpp"

void test_phi_sample()
{
  MyRand myrand;
  
  auto u = myrand.u();
  double A;
  double phis, As;
  sample_phi(u, M_PI_4, phis, A, As);
  
  for (double phi = M_PI_4; phi >= 0; phi -= 0.001)
    sample_phi(1, phi, phis, A, As);
  
  double phi_start = M_PI_4;
  double phi = phi_start;
  for (int i = 0; i < 64; ++i)
  {
    phi *= 0.5;
    sample_phi(1, phi,             phis, A, As);
    sample_phi(1, phi_start - phi, phis, A, As);
  }
  
}

void test_r_sample()
{
  const double r_s_min = 1.0;
  for (double phi = 0; phi <= M_PI_4; phi += 0.01)
  {
    std::cout << "\n";
    const double epsm = 1.*std::numeric_limits<double>::epsilon();
    double r_prior = 1. - epsm;
    const double r_s_max = chock_r(phi);
    for (double u = 0; u <= 1; u += 0.01)
    {
      double r_s;
      sample_r(u, phi, r_s);
      
      std::cout << "phi " << phi << ", u " << u << ", r_s " << r_s << std::endl;
      
      assert(r_s >= r_s_min);
      assert(r_s <= r_s_max);
      assert(r_s >= r_prior); // for fixed phi, r_s is increasing by increasing u
      
      
      r_prior = r_s;
    }
  }
}


void sample_chock(int n, double phi, std::vector<double> &phis, std::vector<double> &rs)
{
  phis.resize(n);
  rs.resize(n);
  
  // to get a different distribution each time
  std::random_device r;
  MyRand myrand( r() );

  // to get the same distribution each time.
  // unsigned seed = 12734824;
  // MyRand myrand(seed);
  
  while( n-- > 0)
  {
    auto u_phi = myrand.u();
    auto u_r = myrand.u();
    
    double phi_s, r_s;
    sample_chock(u_phi, u_r, phi, phi_s, r_s);
    
    phis[n] = phi_s;
    rs[n] = r_s;
  }
}


void test_uniform(double phi, int num_samples, std::string fname)
{
  // sample uniformly from some chocks and draw the sample density in a ps file
  // const double phi = M_PI_4; // M_PI/4.;  * 0.25, *0.16667, 0.08
  std::vector<double> phis, rs;
  sample_chock(num_samples, phi, phis, rs);
  
  // convert to cartesian
  const auto n = phis.size();
  std::vector<double> x(n), y(n);
  
  // some test cases
  
  // horizontal
  Point c(0,1), t(0,0), q(tan(phi),0);
  
  // vertical
  // Point c(0,0), t(1,0), q(1,tan(phi));
  
  // rc = radius of Poisson-disk
  // c = center of Poisson-disk
  // is q cw or ccw from t
  
  // double rc = 0.5;
  // Point c(0.5,0.5);
  
  // double rc = 2.0;
  // Point c(-1.3,-0.5);
  
  // double rc = 0.9;
  // Point c(0.1,0.2);
  
  // bool cw = false; // try both
  
  // 0, M_PI/4, 3*M_PI/4, M_PI, 5*M_PI/4., 6*M_PI/4.; 7*M_PI/4., M_PI/6
  //  double t_angle = M_PI * 0.0;
  //  double t_angle = M_PI * 0.2;
  //  double t_angle = M_PI * 0.4;
  //  double t_angle = M_PI * 0.6;
  //    double t_angle = M_PI * 0.8;
  //  double t_angle = M_PI * 1.0;
  //  double t_angle = M_PI * 1.2;
  //  double t_angle = M_PI * 1.4;
  //  double t_angle = M_PI * 1.6;
  //  double t_angle = M_PI * 1.8;
  //  double t_angle = M_PI * 2.0;
  
  //  Point t( rc*cos(t_angle), rc*sin(t_angle) );
  //  double h = tan(phi);
  //  Point q(0,0);
  //  if (cw)
  //  {
  //    q.x = t.x - h*t.y;
  //    q.y = t.y + h*t.x;
  //  }
  //  else
  //  {
  //    q.x = t.x + h*t.y;
  //    q.y = t.y - h*t.x;
  //  }
  //
  //  t+=c;
  //
  //  q+=c;
  //
  
  
  
  for (size_t i=0; i<n; ++i)
  {
    chock_to_cartesian(c,t,q, phis[i], rs[i], x[i], y[i]);
  }
  
  fname += ".ps";
  std::cout << "writing uniform chock sampling to " << fname << std::endl;
  std::ofstream fout;
  fout.open (fname);
  
  // display in ps file
  fout << "%!PS\n"
  << "0.035 0.000 0.000 setrgbcolor\n"
  //  << "72 72 scale\n"
  //  << "4.25 5.5 translate\n"
  << "144 144 scale\n"
  << "2.125 2.75 translate\n"
  << "0.0005 setlinewidth\n";
  
  // bounding box
  fout << "\n% bounding box\n";
  fout << "newpath\n" <<
  "-1.5 -1.5 moveto\n" <<
  " 1.5 -1.5 lineto\n" <<
  " 1.5  1.5 lineto\n" <<
  " -1.5  1.5 lineto\n" <<
  "closepath\n" << "stroke\n";
  
  plot_chock(c,t,q, fout);
  
  // draw the Points in the chock
  fout << "% sample points in chock\n";
  for (size_t i=0; i<n; ++i)
  {
    fout << "newpath ";
    fout << x[i] << " " << y[i] <<
    " 0.001 0 360 arc closepath fill\n";    
  }
  
  fout << "\nshowpage\n";
  fout.close();
}

void chock_sample_figure()
{
  test_uniform(M_PI_4, 4000, "sc1");
}

