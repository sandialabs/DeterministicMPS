//
//  TestSimpleMPS.cpp
//  AvoidReject
//
//  Created by Mitchell, Scott A on 12/8/21.
//  Copyright © 2021 Mitchell, Scott A. All rights reserved.
//

#include <iostream>
#include <algorithm>
#include <ctime>
#include <chrono>

#include "TestSimpleMPS.hpp"
#include "SimpleMPS.hpp"

void MPStime()
{
  std::cout <<"\nMPStime()\n";
  SimpleMPS mps;
  mps.make_grid(600, 900);
  mps.periodic = true;
  // mps.periodic = false;
  mps.do_progress = false;
  mps.do_plot = false;
  mps.do_extra_checks = false;

  mps.sweep_grid();
}

void MPStimelong()
{
  std::cout <<"\nMPStimelong()\n";
  SimpleMPS mps;
  mps.make_grid(2.25*600, 2.25*900);
  mps.periodic = true;
  // mps.periodic = false;
  mps.do_progress = false;
  mps.do_plot = false;
  mps.do_extra_checks = false;

  mps.sweep_grid();
}

void MPS0()
{
  std::cout <<"\nMPS0()\n";
  SimpleMPS mps;

  // mps.make_grid(6, 9);
  mps.make_grid(60, 90);
  mps.periodic = true;
  // mps.periodic = false;
  mps.do_progress = true;
  mps.do_plot = true;
  mps.do_extra_checks = true;

  auto numsamples = mps.sweep_grid();
  std::cout << "SimpleMPS numsamples = " << numsamples << std::endl;
  mps.plot_grid("MPS", 0);
  mps.plot_dots("MPS", 0);
}

// verify output via grid_covered and samples_separated
// do a timing study for comparison


void MPSscaling(bool periodic, int maxt)
{

  std::cout <<"\nMPSscaling() linear time comparison";
  std::cout <<" " << (periodic? "periodic" : "chopped") << "\n";
  std::vector<int> am = {8, 12, 18, 27, 40, 60, 90, 135, 202, 304,  456, 684, 1026, 1539, 2308, 3462, 5193, 7790, 11684}; //19
  std::vector<int> an = {8, 12, 18, 27, 41, 60, 90, 135, 203, 303,  456, 684, 1026, 1539, 2309, 3463, 5193, 7789, 11684};
  std::vector<long> ag;
  std::cout << "Series of grid cells: ";
  for (int i = 0; i< (int) am.size(); ++i)
  {
    ag.push_back( (long) am[i] * (long) an[i] );
    std::cout << ag.back() << " ";
  }
  std::cout << std::endl;
  
  // int maxt = 17; // 8, 12, 16, 20, with 18, the last one takes 5 minutes and makes 21M samples, slows down due to memory issues?
  maxt = std::min( {maxt, (int) am.size(), (int) an.size()} );

  std::cout << "running test sizes 0..." << maxt-1 << std::endl;
  std::cout << " t   cells  samples  cpu_seconds  cell/s samp/s  cell/samp" << std::endl;

  for (int t=0; t<maxt && t<(int)am.size() && t<(int)an.size(); ++t)
  {
    // use std::clock() on xcode and the mac
    std::clock_t startcputime = std::clock();
    // use steady_clock on linux for greater accuracy.
    // auto startcputime = std::chrono::steady_clock::now();
    
    int gridm = am[t];
    int gridn = an[t];

    SimpleMPS grid;
    unsigned seed = 2453417324;
    grid.myrand = MyRand(seed);
    grid.periodic = periodic;
    grid.do_extra_checks = false;
    
    grid.make_grid(gridm, gridn);
    auto num_samples = grid.sweep_grid();
    
    double cpu_duration = (std::clock() - startcputime) / (double)CLOCKS_PER_SEC;
    // auto endcputime = std::chrono::steady_clock::now();
    // double cpu_duration = std::chrono::duration_cast<std::chrono::duration<double> >(endcputime - startcputime).count();
    
    // printf("%d Cells = %ld. Cells per second processed = %g\n", t, ag[t], ag[t] / cpu_duration );
    // std::cout.flush();
    
    printf("%2d %8ld %8ld   %10.3e %7.0f %6.0f  %5.3f\n", t, ag[t], num_samples, cpu_duration, ag[t] / cpu_duration, num_samples / cpu_duration, (double)ag[t]/ (double)num_samples );
    // std::streamsize ss = std::cout.precision(4);
    //std::cout << ag[t] << " \t" << num_samples << " \t" << cpu_duration << " \t" << ag[t] / cpu_duration << " \t" << num_samples / cpu_duration << std::endl;
    //std::cout.precision(ss);
  }
  std::cout << std::endl;
}


void SimpleMPS_runtime_table()
{
  MPSscaling(true, 17);
}

void SimpleMPS_sampling(int gridm, int gridn, bool periodic, unsigned int rand_seed_p)
{
  std::cout << "SimpleMPS_sampling\n";
  std::clock_t startcputime = std::clock();

  SimpleMPS grid;
  grid.myrand = MyRand(rand_seed_p);
  grid.periodic = periodic;
  grid.do_extra_checks = false; // = true;
  grid.do_progress = true; 

  grid.make_grid(gridm, gridn);
  auto num_samples = grid.sweep_grid();
  
  double cpu_duration = (std::clock() - startcputime) / (double)CLOCKS_PER_SEC;
  long num_cells = gridm*gridn;
  
  std::cout << "cells  samples  cpu_seconds  cell/s samp/s  cell/samp" << std::endl;
  printf("%8ld %8ld   %10.3e %7.0f %6.0f  %5.3f\n", num_cells, num_samples, cpu_duration, num_cells / cpu_duration, num_samples / cpu_duration, (double)num_cells/ (double)num_samples );
}

void SimpleMPS_stats()
{
  int n=10, m=10;
  while (n*m < 24000000)
  {
    SimpleMPS_sampling(n,m);
    n *= sqrt(2.);
    m =n;
  }
}

