//
//  Grid.hpp
//  AvoidReject
//
//  Created by Mitchell, Scott A on 11/12/21.
//  Copyright © 2021 Mitchell, Scott A. All rights reserved.
//

#ifndef Grid_hpp
#define Grid_hpp

#include <vector>
#include <cassert>
#include "Geometry.hpp"
#include "Triangulate.hpp" // for class Triangles
#include "MyRand.hpp"

#include <set>

class GridSquare
{
public:
  // i,j location in Grid is not stored
  //  int i, j;
  //  double x0() const {return (double)i*sqrt(0.5);}
  //  double y0() const {return (double)j*sqrt(0.5);}
  //  double x1() const {return (double)(i+1)*sqrt(0.5);}
  //  double y1() const {return (double)(j+1)*sqrt(0.5);}

  // Poisson-point in square.
  //   p is the point in the global coordinate system 
  Point p;
  double time = 0.;
  int num_earlier = 0;

  // true means there is no place left for a sample. If it has a sample, it is covered.
  bool is_covered = false;

  // true if we've accepted a sample in this square
  bool has_sample = false;

  // is_early means it's been identified as locally early and is in the queue
  bool is_early = false;
  
  // true means we haven't built the poly for this square yet.
  bool poly_uninit = true;

  std::vector<Polygon> polys;

  // disks that might intersect the square, but don't cover the candidate sample so far
  std::vector<Point> deferred_disks;
  
  // set the flags that this square's sample is accepted
  void set_accepted();
  
  // StatusFn of point in square, exactly one should be true
  bool accepted() const {return has_sample;}
  bool early() const {return is_early && !has_sample;}
  bool empty() const {return !has_sample && is_covered;}
  bool not_early() const {return !is_early && !has_sample && !is_covered;} 

  // debug
  // return true if p lies inside the square's polygons
  // i.e. false if p is stricly outside the square, or inside a circle
  bool check_point_in_poly(int i, int j) const;

  // debug
  std::set<int> my_earlier_squares;
  // int times_resampled = 0;
};
typedef bool (GridSquare::*StatusFn)(void) const;

// common stuff for a grid, may be used by SimpleMPS, etc.
// true if p is inside the polys
//   the exact test is just that it is outside all disks bounding the poly and inside the square i,j
bool check_point_in_poly(int i, int j, const Point &p, const std::vector<Polygon> &polys);
Point random_point_in_square(int i, int j, MyRand &myrand);
bool disk_covers_square( const Point &c, int i, int j );
bool disk_intersects_square( const Point &c, int i, int j );

// grid navigation
void neighbors(const int gridm, const int gridn, bool periodic, const int k, 
               std::vector<int> &ks, std::vector<Point> *np_to_p = nullptr);
int TwoDtoOneD(int gridm, int gridn, int i2, int j2);
void OneDtoTwo(int gridm, int gridn, int k1, int &i2, int& j2);

// i/o
double compute_scale(int gridm, int gridn, double &dotsize);
void grid_preamble(double scale, std::ofstream &fout);
void plot_grid_boundary(int gridm, int gridn, std::ofstream &fout);
void plot_squares(int gridm, int gridn, std::ofstream &fout);

struct TriangulatedPoly
{
  double poly_area;
  Chocks chocks;
  Loop loop;
  Triangles triangles;
};
typedef std::vector<TriangulatedPoly> SampleWorkspace;
void sample(int i, int j, std::vector<Polygon> &polys, MyRand &myrand, SampleWorkspace &sample_workspace, bool do_simple_triangulation,
            int &which_poly, Point &p, double &area, int &error_code);

class Grid
{
public:
  // interface

  // algorithm parameters
  
  // domain is a torus if true, else a rectangle
  bool periodic = false;
  
  // random
  MyRand myrand;  // seeds uniform distribution for positions
  std::mt19937 gen; // seeds exponential distribution for times

  // simple or complicated triangulation criteria
  //   faster runtime or prettier pictures?
  bool do_simple_triangulation = true;
  
  // number of "prepasses" where we save a little time by not incrementing Grid::num_earlier (antecedents)
  int point_prepasses = 7; // default 7
  int square_prepasses = 0; // default 0
  bool do_point_neighbors = false; // default false
  
  // methods
  
  // Initialize the grid
  //   number of grid squares in x and y: m x-i values, n y-j values
  void make_grid(int m, int n);

  // Do the sampling.
  // Return the number of accepted samples.
  long sweep_grid();

  // query
  int m() const {return gridm;}
  int n() const {return gridn;}

  // i/o control
  bool do_plot = false;
  double plot_dotsize_factor = 1.;
  double plot_linewidth_factor = 1.;bool do_print = false;
  bool   do_plot_sampletime = false;
  bool do_extra_checks = true;
  bool do_progress = false;
  bool do_counts = false;
  int do_PSA_points = 0; // if >=1, will be written to file "pointsK.txt", where K is the value of do_PSA_points
  bool do_early_stats = false;
  bool do_iter_timing = false;
  
  // i/o
  void print_times() const;
  void print_early() const;
  void print_early_count() const;
  void plot_grid(std::string fname, int id) const;
  void plot_dots(std::string fname, int id) const;
  void plot_gridsquare(int i, int j, std::string fname) const;
  void plot_gridsquare(int i, int j, const Chocks &chocks, Loop &loop, const Triangles &triangles, std::string fname) const;
  void print_PSA_points(std::string fname, long num_accepted = -1 );
  void write_early(int iter);

private:
  // data
  int gridm, gridn;
  std::vector<GridSquare> grid_squares;

  // methods
  
  // get square
  GridSquare &grid_square(int i, int j);
  
  // methods exposed for debugging only
public:
  const GridSquare &const_grid_square(int i, int j) const;
  int TwoDtoOneD(int i2, int j2) const;
  void OneDtoTwo(int k1, int &i2, int& j2) const;
private:

  SampleWorkspace sample_workspace;
  void resample(int k);

  // random time
  double random_expovariate(double area);
  
  // scoop disk out of grid_square[k]
  //   return true if the candidate sample of k changed
  ScoopWorkspace scoop_workspace;
  bool scoop_square(int k, Point &disk);

  // fill ks with indices of square-neighbors
  //   optional np_to_p is the coordinate transform from the neighbor square to k
  //     empty unless periodic and close to boundary.
  void neighbors(const int k, std::vector<int> &ks, std::vector<Point> *np_to_p = nullptr) const;

  // point-based earlier methods
  // true if the square k is blocked from accepting k's sample by square[nk]
  bool point_earlier( int k, int nk, int ni, std::vector<Point> &np_to_p, bool flip_sign ) const;
  bool locally_point_early(int k);
  void find_point_early(std::vector<int> &early_squares);
  void count_point_early(std::vector<int> &early_squares); // increments num_earlier
  void update_early_incremental_pointneighbors(int k, Point& old_p, double old_time, bool old_covered, std::vector<int> &potential_early_squares );

  // square-based earlier methods
  // g is earlier than n : bool g_is_first =  g.time < n.time || (g.time == n.time && nk < gk);
  bool locally_early(int k);
  void update_early_incremental_squareneighbors(int k, double old_time, std::vector<int> &early_squares );
  void find_early(std::vector<int> &early_squares);
  void count_early(std::vector<int> &early_squares); // increments num_earlier
  void zero_early(); // resets count to zero, needed after "find_early" and the counts

  // accept sample in square i,j
  struct AcceptWorkspace
  {
    std::vector<int> ks; //neighbor indices
    std::vector<Point> np_to_p;
    std::vector<int> potential_early_squares;
  };
  AcceptWorkspace accept_workspace;
  void accept_sample_update_pointneighbors(int k, std::vector<int> &early_squares);  // with updated num_earlier counts
  void accept_sample(int k, std::vector<int> &early_squares); // with updated num_earlier counts, _update_squareneighbors
  void accept_sample_no_counts(int k); // agnostic to neighbor type
  
  // find early squares and accept their samples
  // return the number of accepted samples
  long prepass_early(std::vector<int> &early_squares);
  long prepass_point_early(std::vector<int> &early_squares, int pass);

  // i/o
  // plot the samples whose status_fn is true in the given rgb_color
  //   return the number of such samples
  int plot_samples_in_color( StatusFn status_fn, std::string rgb_color, std::string rgb_color_periodic,
                            bool do_circle, std::ofstream &fout, double dotsize = 0.02, double linefactor=1.0) const;
  
  // verification
  void check_status() const;
  bool samples_separated(bool do_progress) const;
  bool grid_covered(bool do_progress) const;
  bool closest_other_disk(int k0, int k1, Point &p, double &dist2 ) const;
  
  // void report_resample_stats();
  void plot_samples(std::ofstream &fout, double dotsize) const;

  friend class TestGrid;
};

double compute_scale(int gridm, int gridn, double &dotsize);
void plot_g(const Point &p, double time, const std::vector<Polygon> &polys, double scale, double dotsize,
            int i, int j, std::string fname);
void plot_g(const Point &p, double time, double scale, double dotsize,
            int i, int j, const Chocks &chocks, Loop &loop, const Triangles &triangles, std::string fname);
// return true if p is close to the x border, and return p's ghost in ghost_p
bool ghost_x( int m, int /*n*/, const Point &p, Point &ghost_p );
bool ghost_y( int /*m*/, int n, const Point &p, Point &ghost_p );
bool ghost_xy( int m, int n, const Point &p, Point &ghost_p );
void fout_samp(int i, int j, const Point &p, double time, bool do_circle, bool do_plot_sampletime, double dotsize, double dotfactor, std::string preamble, std::ofstream &fout);
// get indices i,j of the grid square containing p.
//    wrap p to be in the grid as if it were periodic
void square_containing_point(Point &p, int m, int n, bool periodic, int&i , int&j);


// ==== implementations
inline
Point random_point_in_square(int i, int j, MyRand &myrand)
{
  const auto a = sqrt(0.5);
  return Point( (i+myrand.u())*a, (j+myrand.u())*a );
}

inline
void GridSquare::set_accepted()
{
  has_sample = true;
  is_covered = true;
  polys.clear();
  // time = -time; // for convenience debugging, but updating point-neighbors needs the old time
}

inline
int Grid::TwoDtoOneD(int i2, int j2) const
{
  assert( i2>=0 && i2<gridm);
  assert( j2>=0 && j2<gridn);
  return i2 + j2*gridm;
}

inline
int TwoDtoOneD(int gridm, int gridn, int i2, int j2)
{
  assert( i2>=0 && i2<gridm);
  assert( j2>=0 && j2<gridn);
  return i2 + j2*gridm;
}

inline
void Grid::OneDtoTwo(int k1, int &i2, int& j2) const
{
  assert(k1>=0 && k1 < (int) grid_squares.size());
  i2 = k1 % gridm;
  j2 = k1 / gridm;
}

inline
void OneDtoTwo(int gridm, int gridn, int k1, int &i2, int& j2)
{
  assert(k1>=0 && k1 < gridm*gridn);
  i2 = k1 % gridm;
  j2 = k1 / gridm;
}


inline
GridSquare &Grid::grid_square(int i, int j)
{
  const auto k = TwoDtoOneD(i,j);
  assert(k<(int)grid_squares.size());
  return grid_squares[k];
}

inline
const GridSquare &Grid::const_grid_square(int i, int j) const
{
  const auto k = TwoDtoOneD(i,j);
  assert(k<(int)grid_squares.size());
  return grid_squares[k];
}

inline
bool disk_covers_square( const Point &c, int i, int j )
{
  const double a = sqrt(0.5);
  // find distance to farthest corner
  const double dx = (c.x>a*(i+0.5)) ? c.x-a*i : c.x-a*(i+1);
  const double dy = (c.y>a*(j+0.5)) ? c.y-a*j : c.y-a*(j+1);
  return (dx*dx + dy*dy <= 1.);
}

inline
bool disk_intersects_square( const Point &c, int i, int j )
{
  const double a = sqrt(0.5);
  // find distance to closest side or corner
  const double dx = (c.x<a*i ? c.x-a*i : (c.x>a*(i+1) ? c.x - a*(i+1) : 0));
  const double dy = (c.y<a*j ? c.y-a*j : (c.y>a*(j+1) ? c.y - a*(j+1) : 0));
  return (dx*dx + dy*dy <= 1.);
}

#endif /* Grid_hpp */
