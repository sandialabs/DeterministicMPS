//
//  SimpleMPS.hpp
//  AvoidReject
//
//  Created by Mitchell, Scott A on 12/8/21.
//  Copyright © 2021 Mitchell, Scott A. All rights reserved.
//

#ifndef simpleMPS_hpp
#define simpleMPS_hpp

// implementation of SimpleMPS in the same framework as avoid-reject
// SimpleMPS = A Simple Algorithm for Maximal Poisson-Disk Sampling in High Dimensions

#include "Grid.hpp"

// count the number of various types of operations
// comment this for runtime studies
// define DO_OPS 1

struct MPSSquare
{
  std::vector<int> subcells;
  Point p;
  bool has_sample = false;
};

typedef std::pair<int,int> Cell;

class SimpleMPS
{
public:
public:
  // interface
  bool periodic = false;
  MyRand myrand;

  void make_grid(int m, int n);
  long sweep_grid();

  // i/o control
  bool do_plot = false;
  bool do_print = false;
  bool do_extra_checks = true;
  bool do_progress = false;
  bool do_counts = false;
  
  // i/o
  void plot_dots(std::string fname, int id) const;
  void plot_grid(std::string fname, int id) const;
  void plot_pool(std::string fname, int id) const;

private:
  int gridm, gridn;
  int refinement_level = 0;
  unsigned twoR = 1;
  static const int level_limit = 16;
  std::vector< MPSSquare > squares;
  std::vector< Cell > pool, pool2;
  
  void random_point( const Cell &c, Point &p );

  // top level grid cell containing subcell
  Cell parent_cell(Cell subcell) const;

  // coordinates of corners of the cell
  void corners(const Cell &c, Point &p00, Point &p10, Point &p01, Point &p11) const;

  // true if cell c is covered by a disk
  // bool covered_cell(Cell c) const; // unused
  // true if point in parent cell k is covered by a disk
  bool covered_point(long k, std::vector<int> &ks, std::vector<Point> &np_to_p) const;

  void refine();
  long sample_grid();

  int plot_samples(std::ofstream &fout, double dotsize, bool do_circle) const;

};


#endif /* simpleMPS_hpp */
