//
//  TimelessMPS.hpp
//  AvoidReject
//
//  Created by Mitchell, Scott A on 1/24/22.
//  Copyright © 2022 Mitchell, Scott A. All rights reserved.
//

#ifndef TimelessMPS_hpp
#define TimelessMPS_hpp

#include "Geometry.hpp"
#include "Grid.hpp"

class TimelessSquare
{
public:
  // Poisson-point in square.
  //   p is the point in the global coordinate system
  Point p;

  // true means there is no place left for a sample. If it has a sample, it is covered.
  bool is_covered = false;

  // true if we've accepted a sample in this square
  bool has_sample = false;

  std::vector<Polygon> polys;
  double area = 0.; // 0.5;

  // set the flags that this square's sample is accepted
  void set_accepted();
  
  // StatusFn of point in square, exactly one should be true
  bool accepted() const {return has_sample;}
  bool empty() const {return !has_sample && is_covered;}

  // debug
  // return true if p lies inside the square's polygons
  // i.e. false if p is stricly outside the square, or inside a circle
  bool check_point_in_poly(int i, int j) const;
};

class TimelessGrid
{
public:
  // interface

  // algorithm parameters
  
  // domain is a torus if true, else a rectangle
  bool periodic = false;
  
  // random
  MyRand myrand;  // seeds uniform distribution for positions

  // simple or complicated triangulation criteria
  //   faster runtime or prettier pictures?
  bool do_simple_triangulation = true;
  
  // methods
  
  // Initialize the grid
  //   number of grid squares in x and y: m x-i values, n y-j values
  void make(int m, int n);

  // Do the sampling.
  // Return the number of accepted samples.
  long sweep();

  // query
  int m() const {return gridm;}
  int n() const {return gridn;}

  // i/o control
  bool do_plot = false;
  double plot_dotsize_factor = 1.;
  double plot_linewidth_factor = 1.;bool do_print = false;
  bool do_extra_checks = true;
  bool do_progress = false;
  bool do_counts = false;
  int do_PSA_points = 0; // if >=1, will be written to file "pointsK.txt", where K is the value of do_PSA_points
  bool do_iter_timing = false;
  
  // i/o
  void plot_grid(std::string fname, int id) const;
  void plot_dots(std::string fname, int id) const;
  void plot_square(int i, int j, const Chocks &chocks, Loop &loop, const Triangles &triangles, std::string fname) const;
  void plot_square(int i, int j, std::string fname) const;
  //  void print_PSA_points(std::string fname, long num_accepted = -1 );
  
private:
  // data
  int gridm, gridn;
  std::vector<TimelessSquare> squares;

  // tree for selecting next square uniform by area
  std::vector< std::vector<double> > tree;
  void initialize_tree();
  void update_area(int k, double new_area);
  int random_square();
  // return the square corresponding to the given area (some fraction of the entire tree)
  int treewalk(double area);

  
  // methods
  
  
  // methods exposed for debugging and setting up test problems only
public:
  TimelessSquare &square(int i, int j);
  const TimelessSquare &const_square(int i, int j) const;
  const TimelessSquare &const_square(int k) const;
  int TwoDtoOneD(int i2, int j2) const;
  void OneDtoTwo(int k1, int &i2, int& j2) const;
private:

  SampleWorkspace sample_workspace;
  // return area of void
  double resample(int k);
  
  // accept sample in square i,j
  struct AcceptWorkspace
  {
    std::vector<int> ks; //neighbor indices
    std::vector<Point> np_to_p;
  };
  AcceptWorkspace accept_workspace;
  void accept(int k);
  
  // scoop disk out of grid_square[k]
  //   return true if the square is changed, and the new uncovered area of the square
  ScoopWorkspace scoop_workspace;
  bool scoop_square(int k, Point &disk, double &new_area);

  // fill ks with indices of square-neighbors
  //   optional np_to_p is the coordinate transform from the neighbor square to k
  //     empty unless periodic and close to boundary.
  // void neighbors(const int k, std::vector<int> &ks, std::vector<Point> *np_to_p = nullptr) const;
  
  void plot_samples(std::ofstream &fout, double dotsize, bool accepted_dots_only) const;
  
  // verification
  void check_status() const;
  
};


// ==== implementations

inline
void TimelessSquare::set_accepted()
{
  has_sample = true;
  is_covered = true;
  polys.clear();
}

inline
int TimelessGrid::TwoDtoOneD(int i2, int j2) const
{
  assert( i2>=0 && i2<gridm);
  assert( j2>=0 && j2<gridn);
  return i2 + j2*gridm;
}

inline
void TimelessGrid::OneDtoTwo(int k1, int &i2, int& j2) const
{
  assert(k1>=0 && k1 < (int) squares.size());
  i2 = k1 % gridm;
  j2 = k1 / gridm;
}

inline
TimelessSquare &TimelessGrid::square(int i, int j)
{
  int k = TwoDtoOneD(i,j);
  return squares[k];
}
inline
const TimelessSquare &TimelessGrid::const_square(int i, int j) const
{
  int k = TwoDtoOneD(i,j);
  return squares[k];
}
inline
const TimelessSquare &TimelessGrid::const_square(int k) const
{
  return squares[k];
}

#endif /* TimelessMPS_hpp */
