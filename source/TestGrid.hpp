//
//  TestGrid.hpp
//  AvoidReject
//
//  Created by Mitchell, Scott A on 11/12/21.
//  Copyright © 2021 Mitchell, Scott A. All rights reserved.
//

#ifndef TestGrid_hpp
#define TestGrid_hpp

#include "Grid.hpp"


class TestGrid
{
public:
  
  static bool do_simple_triangulation; // if false, does the slower better looking triangulation
  static bool do_early_stats;

  // just first steps
  static void test_grid_0(); // does both periodic and bounded
  static void test_grid_1(bool periodic);
  static void test_grid_2(bool periodic);
  
  // full mps
  static void test_grid_3(bool periodic);
  static void test_grid_4(bool periodic);
  static void test_grid_5(bool periodic);
  static void test_grid_6(bool periodic);
  static void test_grid_7(bool periodic);

  // huge number of random grid sizes.
  // 3,653,086,500 cells total.
  // Debug with extra checks on, etc, takes 3 days to run: 1 billion a day, or 5M / hour, 84k/minute, 1,400/second
  // Optimized code, about 30x faster
  static void test_grid_8(bool periodic);

  // demonstrate linear timing
  // do tests 0...maxt
  static void test_grid_9(bool periodic, int maxt = 12);

  // weird triangulation case
  static void test_poly_degeneracy_1();

  // rerun any failures here
  static void test_grid_8_isolate1();
  
  static void test_grid(int gridm, int gridn, bool periodic=true, unsigned int rand_seed_p=345038235, unsigned int rand_seed_time=93845671);

  static void test_grid_early(bool periodic);

private:
  static void print_all_neighbors(const Grid &grid);
  static void check_all_neighbors(const Grid &grid);
  static void print_times(Grid &grid);
  static void print_early(Grid &grid);

  // void check_status(const Grid &grid);
  
};

// figure of a run for the paper
void grid_figure();

// figure of a fourier spectrum for the paper
void PSA_point_set_analysis_figure();

// timing as the problem size increases
void DeterministicMPS_runtime_table();


#endif /* TestGrid_hpp */
