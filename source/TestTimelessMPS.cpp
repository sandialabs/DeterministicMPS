//
//  TestTimelessMPS.cpp
//  AvoidReject
//
//  Created by Mitchell, Scott A on 1/24/22.
//  Copyright © 2022 Mitchell, Scott A. All rights reserved.
//

#include "TestTimelessMPS.hpp"

#include <iostream>
#include <ctime>
#include <chrono>
#include <algorithm>

#include "TimelessMPS.hpp"

void timeless_test_0()
{
  std::cout <<"\ntimeless_test_0() multitest";
  TimelessGrid grid;
  grid.do_plot = true;
  grid.do_progress = true;
  // grid.periodic = false;
  grid.periodic = true;
  grid.make(6,8);
  grid.sweep();
}


void timeless_test_3(bool periodic)
{
  std::cout <<"\ntimeless_test_3()";
  std::cout <<" " << (periodic? "periodic" : "chopped") << "\n";
  TimelessGrid grid;
  int gridm = 6;
  int gridn = 8;
  grid.myrand = MyRand(123409234);
  grid.periodic = periodic;

//  grid.do_print = true;
//  grid.do_counts = true;
//  grid.do_progress = true;
//  grid.do_plot = true;
  grid.do_extra_checks = true;
  // grid.do_plot_sampletime = true;
  
  grid.make(gridm, gridn);
  
  // manually assign times and positions
  const auto a = sqrt(2.) / 2.;
  const auto a2 = a / 2.;
  for (int i=0; i<gridm; ++i)
  {
    for (int j=0; j<gridn; ++j)
    {
      grid.square(i,j).p = Point(i*a+a2,j*a+a2);
    }
  }
  
  std::cout << "corner-sweep manual times" << std::endl;
  
  grid.sweep();
}

void timeless_test_4(bool periodic)
{
  std::cout <<"\ntimeless_test_4()";
  std::cout <<" " << (periodic? "periodic" : "chopped") << "\n";
  TimelessGrid grid;
  int gridm = 6;
  int gridn = 6;
  grid.myrand = MyRand(123409234);
  grid.periodic = periodic;

//  grid.do_progress = true;
//  grid.do_counts = true;
//  // grid.do_iter_timing = true;
//  grid.do_plot = true;
//  grid.do_plot_sampletime = true;
//  grid.do_print = true;
//  grid.do_extra_checks = true;
  
  grid.make(gridm, gridn);
  
  grid.sweep();
}

void timeless_test_5(bool periodic)
{
  std::cout <<"\ntimeless_test_5()";
  std::cout <<" " << (periodic? "periodic" : "chopped") << "\n";
  TimelessGrid grid;
  int gridm = 60;
  int gridn = 90;
  grid.myrand = MyRand(134539234);
  grid.periodic = periodic;

  grid.make(gridm, gridn);
  
  grid.sweep();
}

void timeless_test_6(bool periodic)
{
  std::cout <<"\ntimeless_test_6()";
  std::cout <<" " << (periodic? "periodic" : "chopped") << "\n";
  TimelessGrid grid;
  int gridm = 600;
  int gridn = 900;
  grid.myrand = MyRand(34539234);
  grid.periodic = periodic;
  grid.do_extra_checks = false;

  grid.make(gridm, gridn);
  
  grid.sweep();
}

void timeless_test_7(bool periodic)
{
  std::cout <<"\ntimeless_test_7() 8.64M grid squares";
  std::cout <<" " << (periodic? "periodic" : "chopped") << "\n";
  TimelessGrid grid;
  int gridm = 600 * 4;
  int gridn = 900 * 4;
  grid.myrand = MyRand(24453417);
  grid.periodic = periodic;

  grid.make(gridm, gridn);
  
  grid.sweep();
}

void timeless_test_8(bool periodic)
{
  std::cout <<"\ntimeless_test_8() multitest";
  std::cout <<" " << (periodic? "periodic" : "chopped") << "\n";
  int maxgridm = 600;
  int maxgridn = 900;
  std::cout << " up to (" << maxgridm << "," << maxgridn << ")";
  long numcells = 0;
  for (int gridm = 3; gridm<maxgridm; gridm+=4)
  {
    for (int gridn = 2; gridn<maxgridn; gridn+=5)
    {
      numcells += gridm * gridn;
    }
  }
  std::cout << ", total number of cells = " << numcells;
  for (int gridm = 3; gridm<maxgridm; gridm+=4)
  {
    std::cout << std::endl;
    for (int gridn = 2; gridn<maxgridn; gridn+=5)
    {
      std::cout <<"(" << gridm << "," << gridn << ") ";
      const auto a =sqrt(0.5);
      if (periodic && (gridn*a < 2 || gridm*a < 2))
      {
        // code doesn't handle the case of tiny periodicity where a d/  isk overlaps with itself
        std::cout << "=too small ";
        continue;
      }
      TimelessGrid grid;
      unsigned seed = 24453417+ 613*gridm + 227*gridn;
      grid.myrand = MyRand(seed);
      grid.periodic = periodic;

      grid.make(gridm, gridn);
      grid.sweep();
    }
  }
  std::cout << std::endl;
}

void timeless_test_9(bool periodic, int maxt)
{

  std::cout <<"\ntimeless_test_9() scaling table";
  std::cout <<" " << (periodic? "periodic" : "chopped") << "\n";
  std::vector<int> am = {8, 12, 18, 27, 40, 60, 90, 135, 202, 304,  456, 684, 1026, 1539, 2308, 3462, 5193}; // 17, last one runs out of memory on Scott's 8Gb laptop
  std::vector<int> an = {8, 12, 18, 27, 41, 60, 90, 135, 203, 303,  456, 684, 1026, 1539, 2309, 3463, 5193};
  std::vector<long> ag;
  std::cout << "Series of grid cells: ";
  for (int i = 0; i< (int) am.size(); ++i)
  {
    ag.push_back( (long) am[i] * (long) an[i] );
    std::cout << ag.back() << " ";
  }
  std::cout << std::endl;

  maxt = std::min( {maxt, (int) am.size(), (int) an.size()} );
  
  std::cout << "Beware the runtime tends to bounce around a lot when rerunning: second and subsequent runs tend to be faster.\n"
  "Fixed random number generator seeds are used, so that is not the cause of this variation.\n"
  "The speculation is this is caused by memory allocation, cleanup, the state of the cache, and the state of other processes running on the machine.\n"
  "You can also get different runtimes if you run just a few of the tests at a time :(\n"
  "E.g. Running tests 9 and 10, gives 33k samples/s the first time, then 41, then 49k! This is a 50% increase in efficiency!"  << std::endl;
  std::cout << "running test sizes 0..." << maxt-1 << std::endl;
  std::cout << " t    cells  samples  cpu_seconds  cell/s samp/s  cell/samp" << std::endl;

  // for (int t=11; t<=11; ++t)
  for (int t=0; t<maxt; ++t)
  {
    std::clock_t startcputime = std::clock();
    // auto chrono_start = std::chrono::high_resolution_clock::now();
    int gridm = am[t];
    int gridn = an[t];

    TimelessGrid grid;
    unsigned seed = 2453417324;
    // unsigned seed = 590349834;
    grid.myrand = MyRand(seed);
    grid.periodic = periodic;
    grid.do_counts = false;
    grid.do_extra_checks = false;

    // debug
    // grid.do_counts = true;
    // grid.do_extra_checks = true;
    // grid.do_plot = true;
    // grid.do_progress = true;
    // grid.do_print = true;
    
    grid.make(gridm, gridn);
    auto num_samples = grid.sweep();
    
    auto endcputime = std::clock();
    // auto chrono_end = std::chrono::high_resolution_clock::now();
    double cpu_duration = (endcputime - startcputime) / (double)CLOCKS_PER_SEC;
    // double chrono_duration = std::chrono::duration_cast<std::chrono::duration<double> >(chrono_end - chrono_start).count();
    
    printf("%2d %8ld %8ld   %10.3e %7.0f %6.0f  %5.3f -clock\n", t, ag[t], num_samples, cpu_duration, ag[t] / cpu_duration, num_samples / cpu_duration, (double)ag[t]/ (double)num_samples );
    // printf("%2d %8ld %8ld   %10.3e %7.0f %6.0f  %5.3f -chrono\n", t, ag[t], num_samples, chrono_duration, ag[t] / chrono_duration, num_samples / chrono_duration, (double)ag[t]/ (double)num_samples );
  }
  std::cout << std::endl;
}



void TimelessMPS_runtime_table()
{
  timeless_test_9(true, 16); //16
}
