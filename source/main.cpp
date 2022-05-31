//
//  main.cpp
//  AvoidReject
//
//  Created by Mitchell, Scott A on 7/23/21.
//  Copyright © 2021 Mitchell, Scott A. All rights reserved.
//

#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include <assert.h>

// key methods
#include "ChockSample.hpp"
#include "Geometry.hpp"

// tests
#include "TestChockSample.hpp"
#include "TestGeometry.hpp"
#include "TestTriangulate.hpp"
#include "TestPolygonSample.hpp"
#include "TestGrid.hpp"
#include "TestSimpleMPS.hpp"
#include "TestTimelessMPS.hpp"


void test_chock_sampling()
{
  // options
  
  test_phi_sample();
  
  test_r_sample();
  
  test_uniform();
}

void test_geometry()
{
  test_intersect_circle();
  test_intersect_x();
  test_intersect_y();
}

void test_polygon()
{
  test_polygon0();
  test_polygon1();
  test_polygon2();
  test_polygon3();
  // test_polygon4();
}

void   test_scoop_circle_circle()
{
  test_scoop_circle_circle0();
  test_scoop_circle_circle1();
  test_scoop_circle_circle2();
}

void test_shave()
{
  test_shave_0();
  test_shave_1();
  test_shave_2();
  // test_shave_3(); // long time
}

void test_triangulate()
{
   test_angle_measure();
   test_circum();
  test_triangulate_0();
  test_triangulate_1();
  test_triangulate_2();
 // test_triangulate_3(); // long time
   test_grazing_0();
   test_grazing_1();
  test_grazing_2();
}

void test_polygon_sample()
{
  test_polygon_sample_0();
  test_polygon_sample_1();
  test_polygon_sample_2();
  test_polygon_sample_3();
  // test_polygon_sample_4(); // long time
}

void test_grid()
{
  // fast or pretty?
  // TestGrid::do_simple_triangulation = false;
  TestGrid::do_simple_triangulation = true;

  // just first steps
  // non-periodic
//  TestGrid::test_grid_0(false);
//  TestGrid::test_grid_1(false);
//  TestGrid::test_grid_2(false);

   // periodic
//   TestGrid::test_grid_0(true);
//   TestGrid::test_grid_1(true);
//   TestGrid::test_grid_2(true);

  // actual mps
  
  // non-periodic
  TestGrid::test_grid_3(false);
  TestGrid::test_grid_4(false);
  TestGrid::test_grid_5(false);
  TestGrid::test_grid_6(false); // 10 second timing test
  // TestGrid::test_grid_7(false); // 8.4M samples
  
  // periodic
  TestGrid::test_grid_3(true);
  TestGrid::test_grid_4(true);
  TestGrid::test_grid_5(true);
  TestGrid::test_grid_6(true); // 10 second timing test
  // TestGrid::test_grid_7(true); // 8.4M samples
  
  // test_grid_3 was crash, bad geometry intersection, fixed by changing equality to inequality. Now iter 5, only one ear to cut
  // test_grid_6 was crash, bad resample, because of bad triangulation. Fixed by changing triangulation criteria. Takes 1 second.
  // test_grid_3 had another triangulation failure,
  //   but fixed after changes to triangulation to prefer near-zero angled triangles, and fix the circumcenter calculation for degenerate triangles with two coincident points
  
  // some huge random set of tests
  //  TestGrid::test_grid_8(true); // takes about 3 days to finish
  //  TestGrid::test_grid_8(false); // takes about 3 days to finish
  
  // scaling test
  TestGrid::test_grid_9(false);
  TestGrid::test_grid_9(true);
  
  // TestGrid::test_grid_8_isolate1();

}


void test_scoop_degeneracy_3()
{
  // three circles meeting at a point,
  //  various roundoff cases of how they overlap
  //  include one box side
  //  include two box sides at a corner
  
  
  // scoop
  // triangulate
  // check area calculation is close to zero
  // check triangles are not inside out
  // sample
  // check sample is good, not inside a circle
}

void test_scoop_degeneracy_4()
{
  // four circles, TestGrid::test_grid_3 already covers this.
}

void test_scoop_degeneracy_5()
{
  // five circles
}

void test_scoop_degeneracy_6()
{
  // six circles
}

void debug()
{
  // test_scoop_degeneracy_1();
  TestGrid::test_poly_degeneracy_1();
}

void test_mps()
{
  // MPS0();
  // MPStime();
  // MPStimelong();
  MPSscaling(false);
  MPSscaling(true);
}

void test_early_stats()
{
  TestGrid::do_early_stats = true;
  int n=10, m=10;
  // while (n*m < 1200000) // 1.2M cells takes about 7 seconds on scott's mac
  while (n*m < 12000000) // 12M cells takes about 83 seconds on scott's mac
  {
    TestGrid::test_grid(n,m);
    n *= sqrt(2.);
    m =n;
  }
}

void test_timeless()
{
  timeless_test_0();
  
  timeless_test_3(false);
  timeless_test_4(false);
  timeless_test_5(false);
  timeless_test_6(false);
  // timeless_test_7(false); // big
  timeless_test_8(false);
  timeless_test_9(false, 14);

  timeless_test_3(true);
  timeless_test_4(true);
  timeless_test_5(true);
  timeless_test_6(true);
  // timeless_test_7(true); // big
  timeless_test_8(true);
  timeless_test_9(true, 14);

  TimelessMPS_runtime_table();
}

void make_figures()
{
  // make all the figures in the paper
  
  figure_trim(); // in TestTriangulate
  figure_triangles(); // in TestTriangulate
  
  grid_figure();   // in TestGrid
  
  chock_sample_figure(); // in TestChockSample
  
  DeterministicMPS_runtime_table();
  
  SimpleMPS_runtime_table();
  
  PSA_point_set_analysis_figure(); // requires running PSA 3rd party code later to actually make the figure
}



int main(int argc, const char * argv[])
{
  std::cout.precision(std::numeric_limits<double>::max_digits10);

  make_figures();


  // experiments with variant of algorithm that has no explicit arrival time
  // timeless_test_6(true);
  //
  // timeless_test_9(true, 16);
  // test_timeless();
  
  // TestGrid::test_grid_3(false);
  // TestGrid::test_grid_early(true);
  // test_grid();
  // TestGrid::test_grid_6(true);
  // DeterministicMPS_runtime_table();

  // unit-tests
  //  debug();
  //  
  //  test_chock_sampling();
  //  
  //  test_geometry();
  //  
  //  test_polygon();
  //  
  //  test_scoop_circle_circle();
  //  
  //  test_shave();
  //  
  //  test_triangulate();
  //  
  //  test_polygon_sample();
  //  
  //  test_grid();
  //  
  //  test_mps();
    
  // test_early_stats();
  // SimpleMPS_stats();
  // DeterministicMPS_runtime_table();
  
  return 0;
}
