//
//  Grid.cpp
//  AvoidReject
//
//  Created by Mitchell, Scott A on 11/12/21.
//  Copyright © 2021 Mitchell, Scott A. All rights reserved.
//

#include "Grid.hpp"

#include <iostream>
#include <fstream>
#include <ostream>
#include <iomanip>
#include <algorithm>

#include "Triangulate.hpp"
#include "PolygonSample.hpp"

double Grid::random_expovariate(double area)
{
  std::exponential_distribution<> d(area);
  const auto r = d(gen);
  return r;
}

inline
void random_point_in_point(const Loop &loop, Point &p, MyRand &myrand)
{
  p = loop.front();
  // t = (0.5 + 0.01*myrand.u()) * std::numeric_limits<double>::max();
}

inline
void random_point_in_interval(const Loop &loop, Point &p, MyRand &myrand)
{
  const double u = myrand.u();
  p = loop.front()*u + loop.back()*(1.-u);
  // t = (0.2 + 0.01*myrand.u()) * std::numeric_limits<double>::max();
}

void Grid::make_grid(int m, int n)
{
  gridm = m;
  gridn = n;
  grid_squares.resize(m*n);
  
  // create points with time
  for (int k=0; k<(int)grid_squares.size(); ++k)
  {
    int i, j;
    OneDtoTwo(k,i,j);
    
    auto &g = grid_squares[k];
    g.p = random_point_in_square(i,j,myrand);
    
    // square area is 0.5
    g.time = random_expovariate(0.5);
  }
  accept_workspace.ks.reserve(20);
  accept_workspace.np_to_p.reserve(20);
  accept_workspace.potential_early_squares.reserve(200);
  scoop_workspace.reserve(4); // max of 4 voids is possible
  sample_workspace.resize(4); // max of 4 voids is possible
}

// order of increasing distance to center square
void neighbors(const int gridm, const int gridn, bool periodic, const int k,
               std::vector<int> &ks, std::vector<Point> *np_to_p)
{
  // 5 x 5 subgrid with corners removed, and center i,j removed, centered on (i,j)
  if (np_to_p) np_to_p->clear();
  
  // this square
  int i, j;
  OneDtoTwo(gridm,gridn, k,i,j);
  
  // simple for far-interior
  if (i>=2 && i<gridm-2 && j>=2 && j<gridn-2)
  {
    const auto km1 = k-gridm;
    const auto km2 = km1-gridm;
    const auto kp1 = k+gridm;
    const auto kp2 = kp1+gridm;
    // scan line order
    // ks = {km2-1, km2, km2+1,   km1-2, km1-1, km1, km1+1, km1+2,   k-2,k-1, /*k*/ k+1,k+2,  kp1-2, kp1-1, kp1, kp1+1, kp1+2,  kp2-1, kp2, kp2+1};
    
    // close squares first order
    ks = {km1, k+1, kp1, k-1,   km1-1, km1+1, kp1+1, kp1-1,   km2, k+2, kp2, k-2,  km2-1, km2+1, km1+2, kp1+2, kp2+1, kp2-1, kp1-2, km1-2};
    return;
  }

  // scan line order
  // std::vector<int> is = { i-1,i,  i+1,  i-2,i-1,i  ,i+1,i+2,  i-2,i-1,i+1,i+2,  i-2,i-1,i  ,i+1,i+2,   i-1,i  ,i+1};
  // std::vector<int> js = { j-2,j-2,j-2,  j-1,j-1,j-1,j-1,j-1,  j,j,j,j,          j+1,j+1,j+1,j+1,j+1,   j+2,j+2,j+2};

  // close squares first order
  std::vector<int> is = {i, i+1, i, i-1,   i-1, i+1, i+1, i-1,   i, i+2, i, i-2,  i-1, i+1, i+2, i+2, i+1, i-1, i-2, i-2};
  std::vector<int> js = {j-1, j, j+1, j,   j-1, j-1, j+1, j+1,   j-2, j, j+2, j,  j-2, j-2, j-1, j+1, j+2, j+2, j+1, j-1};
  
  ks.clear();
  
  if (periodic)
  {
    // wrap indices
    for (size_t a = 0; a<is.size(); ++a)
    {
      auto ii = is[a];
      auto jj = js[a];
      if (ii<0)
      {
        ii+=gridm;
      }
      else if (ii>=gridm)
      {
        ii-=gridm;
      }
      if (jj < 0)
      {
        jj+=gridn;
      }
      else if (jj>=gridn)
      {
        jj-=gridn;
      }
      assert(ii>=0 && ii<gridm);
      assert(jj>=0 && jj<gridn);
      ks.push_back( TwoDtoOneD(gridm,gridn, ii,jj) );
    }
    // wrap coordinates
    if (np_to_p)
    {
      np_to_p->resize(20); // zeros
      const double a = sqrt(0.5);
      const double dx = gridm*a;
      const double dy = gridn*a;
      for (size_t b = 0; b<is.size(); ++b)
      {
        auto ii = is[b];
        auto jj = js[b];
        auto &npa = (*np_to_p)[b];
        if (ii<0)
        {
          npa.x -= dx;
        }
        else if (ii>=gridm)
        {
          npa.x += dx;
        }
        if (jj<0)
        {
          npa.y -= dy;
        }
        else if (jj>=gridn)
        {
          npa.y += dy;
        }
      }
    }
  }
  else
  {
    for (size_t a=0; a<is.size(); ++a)
    {
      auto ii = is[a];
      auto jj = js[a];
      if (ii<0 || ii>=gridm || jj<0 || jj>=gridn)
      {
        continue;
      }
      ks.push_back( TwoDtoOneD(gridm,gridn,ii,jj) );
    }
  }
}

void Grid::neighbors(const int k, std::vector<int> &ks, std::vector<Point> *np_to_p) const
{
  return ::neighbors( gridm, gridn, periodic, k, ks, np_to_p);
}

void triangulate_loop(Loop &loop, Triangles &triangles, bool do_simple_triangulation, int &error_code)
{
  if (loop.size()>=3)
  {
    if (do_simple_triangulation)
    {
      triangulate_loop_simple(loop,triangles);
    }
    else
    {
      auto loc_error_code = triangulate_loop(loop,triangles);
      if (loc_error_code)
      {
        loc_error_code = triangulate_loop_simple(loop,triangles);
        if (!loc_error_code)
        {
          std::cout << "Recovered from fancy triangulation warning/error. Retriangulating with a simpler algorithm worked." << std::endl;
        }
        else
        {
          error_code = loc_error_code;
        }
      }
    }
  }
  else
  {
    triangles.clear();
  }
}

void sample(int i, int j, std::vector<Polygon> &polys, MyRand &myrand, SampleWorkspace &sample_workspace, bool do_simple_triangulation,
            int &which_poly, Point &p, double &area, int &error_code)
{
  error_code = 0;

  // pristine square?
  if (polys.empty())
  {
    p = random_point_in_square(i, j, myrand);
    area = 0.5;
    which_poly = -1;
    return;
  }
  // one poly? usualy case, 90%+
  else if (polys.size()==1)
  {
    which_poly = 0;
    const auto min_size = std::min(4, (int) polys.size());
    if (sample_workspace.size() < min_size) sample_workspace.resize(min_size);

    Polygon &poly = polys.front();
    TriangulatedPoly &tp = sample_workspace.front();
    auto &chocks = tp.chocks;
    auto &loop = tp.loop;
    auto &triangles = tp.triangles;
    
    // plot(poly, "square_" + std::to_string(gi) + "_" + std::to_string(gj)); // debug
    shave_chocks(poly, chocks, loop);
    simplify_loop(loop);
    triangulate_loop(loop,triangles,do_simple_triangulation,error_code);
     
    //area
    area = polygon_area(chocks, loop, triangles);
      
    // position
    sample_polygon(chocks, loop, triangles, p, myrand, area);
  }
  // multiple polys
  else
  {
    const auto min_size = std::min(4, (int) polys.size());
    if (sample_workspace.size() < min_size) sample_workspace.resize(min_size);

    // gather poly areas
    area = 0.;
    for (size_t ip = 0; ip<polys.size(); ++ip)
    {
      auto &poly = polys[ip];
      auto &tp = sample_workspace[ip];
      auto &chocks = tp.chocks;
      auto &loop = tp.loop;
      auto &triangles = tp.triangles;
      auto &poly_area = tp.poly_area;
      
      shave_chocks(poly, chocks, loop);
      simplify_loop(loop);
      triangulate_loop(loop,triangles,do_simple_triangulation,error_code);
      poly_area = polygon_area(chocks, loop, triangles);
      area += poly_area;
    }
    
    // select which_poly uniformly by area
    double rand_area = area * myrand.u();
    which_poly = 0;
    while (true)
    {
      rand_area -= sample_workspace[which_poly].poly_area;
      if (rand_area<=0. || which_poly+1 >= polys.size())
      {
        break;
      }
      ++which_poly;
    };

    // sample from polys[which_poly]
    auto &tp = sample_workspace[which_poly];
    auto &chocks = tp.chocks;
    auto &loop = tp.loop;
    auto &triangles = tp.triangles;
    auto &poly_area = tp.poly_area;

    // position (use poly_area)
    sample_polygon(chocks, loop, triangles, p, myrand, poly_area);
  }
}

void Grid::resample(int k)
{
  GridSquare &g = grid_squares[k];
  // gs.times_resampled++;
  
  // shouldn't be calling this on an empty square
  assert( !g.is_covered );
  assert( !g.polys.empty() || g.poly_uninit );

  int i, j;
  OneDtoTwo(k, i,j);
  double polys_area;
  int error_code, which_poly;
  sample(i, j, g.polys, myrand, sample_workspace, do_simple_triangulation,
         which_poly, g.p, polys_area, error_code);
  assert(polys_area>=0.);
  g.time += random_expovariate(polys_area);

  if (error_code)
  {
    std::cout << "resample triangulation issue" << std::endl;

    // redo to step through in debugger
    // auto &poly = g.polys[which_poly];
    auto &tp = sample_workspace[which_poly];
    // auto &chocks = tp.chocks;
    auto &loop = tp.loop;
    // auto &triangles = tp.triangles;
    // auto &poly_area = tp.poly_area;

    Triangles triangles2, triangles3;
    int error_code2 = triangulate_loop(loop,triangles2);
    std::cout << "error_code2=" << error_code2 << std::endl;
    int error_code3 = triangulate_loop_simple(loop,triangles3);
    std::cout << "error_code3=" << error_code3 << std::endl;

  }
  if (do_extra_checks)
  {
    if (!g.check_point_in_poly(i, j))
    {
      std::cout << "error, sample point is not in the poly." << std::endl;
      error_code = 1;
    }
  }
  if (error_code)
  {
    // debug
    auto &poly = g.polys[which_poly];
    auto &tp = sample_workspace[which_poly];
    auto &chocks = tp.chocks;
    auto &loop = tp.loop;
    auto &triangles = tp.triangles;
    // auto &poly_area = tp.poly_area;

    // plot the problem square
    print_poly(poly);
    print_loop(loop);
    plot_gridsquare(i,j,"sample_error");
    plot_gridsquare(i,j,chocks,loop,triangles,"sample_error");
    
    Chocks nochocks;
    Triangles notriangles;
    Loop noloop;
    plot_gridsquare(i,j, chocks, loop, notriangles, "sample_error-chocksonly");
    plot_gridsquare(i,j, nochocks, loop, notriangles, "sample_error-looponly");
    plot_gridsquare(i,j, nochocks, loop, triangles, "sample_error-trianglesonly");
  }
}

bool check_point_in_poly(int i, int j, const Point &p, const std::vector<Polygon> &polys)
{

  // debug, double-check that the sample is in the square
  const double a=sqrt(0.5);
  assert(p.x >= i*a);
  assert(p.x <= (i+1)*a);
  assert(p.y >= j*a);
  assert(p.y <= (j+1)*a);
  const auto t = __DBL_EPSILON__; // side_threshold
  if (p.x<i*a-t || p.x>(i+1)*a+t || p.y<j*a-t || p.y>(j+1)*a+t)
  {
    return false;
  }
  // double-check that the sample is in one of the chocks or polygons, i.e. outside all the circles
  const auto almost_one2 = (1. - 2.*__DBL_EPSILON__)*(1. - 2.*__DBL_EPSILON__); // about circle radius
  for (auto &poly : polys)
  {
    for (auto &s : poly)
    {
      if (s.ctype==CIRCLE)
      {
        const Point v(s.dp - p);
        const auto dist2 = normsquared(v);
        if ( dist2 < almost_one2 )
        {
          std::cout << "candidate sample ";
          print_point(p);
          std::cout <<" is at distances " <<
          std::setprecision(std::numeric_limits<double>::max_digits10) << dist2 <<
          " from an accepted sample ";
          print_point(s.dp);
          std::cout << std::endl;
          // assert(dist2 >= almost_one2);
          return false;
        }
      }
    }
  }
  return true;
  
}

bool GridSquare::check_point_in_poly(int i, int j) const
{
  return ::check_point_in_poly(i,j,p,polys);
}

bool Grid::scoop_square(int k, Point &disk)
{
  auto &g = grid_squares[k];
  int i, j;
  OneDtoTwo(k, i, j);

  assert(!g.is_covered); // should have caught this before calling

  // if the neighbors sample point is still good, defer the bite
  if (point_not_inside_circle(g.p,disk))
  {
    // save circle for biting later as needed
    //   biting will never take place if n_square.is_early, since we'll accept the sample as is
    if (disk_intersects_square( disk, i, j ))
    {
      g.deferred_disks.reserve(4);
      g.deferred_disks.push_back( disk );
    }
    return false;
  }

  // sample is covered - need a new one unless the whole square is empty
  
  // build the square poly if this is the first time
  bool check_disk_covers_polygon = true;
  if (g.poly_uninit)
  {
    g.poly_uninit = false;
    if (disk_covers_square(disk, i, j)) // test timing
    {
      g.is_covered = true;
      g.time = std::numeric_limits<double>::max(); // inf
      return true;
    }
    else
    {
      check_disk_covers_polygon = false;
      const auto a = sqrt(0.5);
      // g.polys.push_back( square_to_poly( (double) i*a, (double) j*a ) );
      g.polys.reserve(4);
      g.polys.emplace_back( square_to_poly( (double) i*a, (double) j*a ) );
    }
  }

  // scoop
  scoop_polys(g.polys, disk, scoop_workspace, check_disk_covers_polygon);

  // scoop deferred disks
  if (!g.deferred_disks.empty())
  {
    for (auto &d : g.deferred_disks)
    {
      scoop_polys(g.polys, d, scoop_workspace, true);
    }
    g.deferred_disks.clear();
  }

  // resample if not covered
  if (g.polys.empty())
  {
    g.is_covered = true;
    g.time = std::numeric_limits<double>::max(); // inf
  }
  else
  {
    resample(k);
  }
  return true;
}

void Grid::accept_sample_no_counts(int k)
{
  bool do_debug = false;
  if (/* DISABLES CODE */ (0) && k==13)
  {
    plot_grid("accepting_sample_", 13 );
    plot_dots("accepting_sample_", 13 );
    do_debug = true;
  }

  // update grid square k
  auto &g = grid_squares[k];
  g.set_accepted();

  // accepted disk, updated later as needed to be in coordinate system of neighbor
  Point disk(g.p);

  //== update neighbors

  // gather neighbors
  //   optimization: re-use storage, move to class so we don't keep making all these vectors
  std::vector<int> &ks = accept_workspace.ks;
  std::vector<Point> &np_to_p = accept_workspace.np_to_p;
  neighbors(k, ks, &np_to_p);
  
  // bite neighborsneighbors
  for (int ki=0; ki< (int) ks.size(); ++ki)
  {
    auto nk = ks[ki];
    auto &ng = grid_squares[nk];
    // skip squares covered previously, e.g. has_sample, or empty polys
    if (ng.is_covered)
    {
      continue;
    }

    if (do_debug)
    {
      int ni, nj;
      OneDtoTwo(nk,ni,nj);
      plot_gridsquare(ni, nj, "neighbor_13_" + std::to_string(nk) + " (" + std::to_string(ni) + "," + std::to_string(nj) + ")" );
    }
    
    // translate accepted disk in coordinate system of neighbor
    if (!np_to_p.empty()) disk = (g.p - np_to_p[ki]);
    scoop_square(nk, disk);
  }
}


void Grid::accept_sample(int k, std::vector<int> &early_squares) // _update_squareneighbors
{
  // update grid square k
  auto &g = grid_squares[k];
  g.set_accepted();

  // accepted disk, updated later as needed to be in coordinate system of neighbor
  Point disk(g.p);

  //== update neighbors

  // gather neighbors
  //   optimization: re-use storage, move to class so we don't keep making all these vectors
  std::vector<int> &ks = accept_workspace.ks;
  std::vector<Point> &np_to_p = accept_workspace.np_to_p;
  neighbors(k, ks, &np_to_p);
  
  // bite neighbors
  for (int ki=0; ki< (int) ks.size(); ++ki)
  {
    auto nk = ks[ki];
    auto &ng = grid_squares[nk];
    if (ng.is_covered)
    {
      continue;
    }

    // tell neighbor square that g is no longer a blocker, nothing to check
    ng.num_earlier--;
    if (do_extra_checks)
    {
      assert( ng.my_earlier_squares.find(k) != ng.my_earlier_squares.end() );
      ng.my_earlier_squares.erase(k);
    }

    // dt = disk-translate, accepted disk in coordinate system of neighbor
    if (!np_to_p.empty()) disk = (g.p - np_to_p[ki]);

    const auto old_time = ng.time;
    if (scoop_square(nk, disk))
    {
      update_early_incremental_squareneighbors( nk, old_time, early_squares);
    }
    // always check, as g might have been the only blocker
    if (ng.num_earlier==0 && !ng.is_early && !ng.is_covered)
    {
      ng.is_early=true;
      early_squares.push_back(nk);
    }
  }
}

void Grid::accept_sample_update_pointneighbors(int k, std::vector<int> &early_squares)
{
  // update grid square k
  auto &g = grid_squares[k];
  g.set_accepted();

  // accepted disk, updated later as needed to be in coordinate system of neighbor
  Point disk(g.p);

  //== update neighbors

  // gather neighbors
  //   optimization: re-use storage, move to class so we don't keep making all these vectors
  std::vector<int> &ks = accept_workspace.ks;
  std::vector<Point> &np_to_p = accept_workspace.np_to_p;
  neighbors(k, ks, &np_to_p);
  
  auto &potential_early_squares = accept_workspace.potential_early_squares;
  
  // bite neighbors
  for (int ki=0; ki< (int) ks.size(); ++ki)
  {
    auto nk = ks[ki];
    auto &ng = grid_squares[nk];
    if (ng.is_covered)
    {
      continue;
    }

    // old point and time
    Point old_p = ng.p;
    double old_time = ng.time;
    bool old_covered = ng.is_covered;

    // dt = disk-translate, accepted disk in coordinate system of neighbor
    if (!np_to_p.empty()) disk = (g.p - np_to_p[ki]);

    if (scoop_square(nk, disk))
    {
      update_early_incremental_pointneighbors( nk, old_p, old_time, old_covered, potential_early_squares);
    }
    
    if (!ng.is_covered)
    {
      // update count if nk used to be blocked by k
      g.is_covered=false;
      std::swap(ng.p,old_p);
      std::swap(ng.time,old_time);
      std::swap(ng.is_covered,old_covered);
      const bool was_nk_blocked_by_k = point_earlier(nk, k, ki, np_to_p, true);
      g.is_covered=true;
      std::swap(ng.p,old_p);
      std::swap(ng.time,old_time);
      std::swap(ng.is_covered,old_covered);
      if (was_nk_blocked_by_k)
      {
        if (do_extra_checks)
        {
          assert( ng.my_earlier_squares.find(k) != ng.my_earlier_squares.end() );
          ng.my_earlier_squares.erase(k);
        }

        if (--ng.num_earlier==0 && !ng.is_early)
        {
          ng.is_early=true;
          early_squares.push_back(nk);
        }
        assert( ng.num_earlier >= 0 );
      }
    }
  }
  // check for early squares once *all* increments have been done
  for (auto &ki : potential_early_squares)
  {
    auto &nn = grid_squares[ki];
    if (nn.num_earlier==0 && !nn.is_early && !nn.is_covered)
    {
      early_squares.push_back(ki);
      nn.is_early=true;
    }
  }
}

bool Grid::locally_early(int gk)
{
  GridSquare &g = grid_squares[gk];
  if (g.is_covered)
  {
    return false;
  }
  std::vector<int> ks;
  neighbors(gk, ks);
  for (int nk : ks)
  {
    auto &n = grid_squares[nk];
    // not early if a neighbor square has a candidate sample with an earlier time
    //   in case of time tie, square with a smaller k is considered earlier
    if (!n.is_covered && ( n.time < g.time || (n.time == g.time && nk < gk) ) )
    {
      return false;
    }
  }
  return true;
}

bool Grid::locally_point_early(int gk)
{
  GridSquare &g = grid_squares[gk];
  if (g.is_covered)
  {
    return false;
  }
  std::vector<int> ks;
  std::vector<Point> np_to_p;
  neighbors(gk, ks, &np_to_p);
  for (int ki=0; ki<(int)ks.size(); ++ki)
  {
    int nk = ks[ki];
    auto &n = grid_squares[nk];
    if (!n.is_covered && point_earlier(gk, nk, ki, np_to_p, false ))
    {
      return false;
    }
  }
  return true;
}

void Grid::update_early_incremental_squareneighbors(int k, double old_time, std::vector<int> &early_squares )
{
  // save count of neighbors that are earlier and later
  int i,j;
  OneDtoTwo(k, i,j);
  GridSquare &g = grid_squares[k];

  std::vector<int> ks;
  neighbors(k, ks);
  
  for (auto nk: ks)
  {
    auto &n = grid_squares[nk];

    if (n.is_covered)
    {
      continue;
    }
    const bool n_was_early =                   ( n.time < old_time || (n.time == old_time && nk < k) );
    const bool n_now_early = (g.is_covered) || ( n.time <   g.time || (n.time ==   g.time && nk < k) );
    
    if (!n_was_early && n_now_early)
    {
      if (do_extra_checks)
      {
        assert( g.my_earlier_squares.find(nk) == g.my_earlier_squares.end() );
        g.my_earlier_squares.insert(nk);
        assert( n.my_earlier_squares.find(k) != n.my_earlier_squares.end() );
        n.my_earlier_squares.erase(k);
      }
      g.num_earlier++;
      n.num_earlier--;
      if (n.num_earlier==0)
      {
        early_squares.push_back(nk);
        assert( !n.is_early );
        n.is_early = true;
      }
    }
    else if (n_was_early && !n_now_early)
    {
      // since we currently only call this when square g has been resampled, not when ng is resampled, this case is never hit
      assert(0);
      if (do_extra_checks)
      {
        assert( g.my_earlier_squares.find(nk) != g.my_earlier_squares.end() );
        g.my_earlier_squares.erase(nk);
        assert( n.my_earlier_squares.find(nk) == n.my_earlier_squares.end() );
        n.my_earlier_squares.insert(k);
      }
      g.num_earlier--;
      n.num_earlier++;
      assert( !g.is_early );
      assert( !g.is_covered);
      if (g.num_earlier==0)
      {
        early_squares.push_back(k);
        g.is_early = true;
      }
    }
  }
}


void Grid::update_early_incremental_pointneighbors(int k, Point &old_p, double old_time, bool old_covered, std::vector<int> &potential_early_squares )
{
  // increment count of neighbors that are earlier, after changing the point and time of square k
  GridSquare &g = grid_squares[k];

  std::vector<int> ks;
  std::vector<Point> np_to_p;
  neighbors(k, ks, &np_to_p);
  
  for (int ni=0; ni< (int)ks.size(); ni++)
  {
    auto nk = ks[ni];
    auto &ng = grid_squares[nk];
    if (ng.is_covered) // e.g. has sample
    {
      continue;
    }
    
    // prior and current state of who is blocked
    std::swap(g.p,old_p);
    std::swap(g.time,old_time);
    std::swap(g.is_covered,old_covered);
    bool was_k_blocked_by_n = point_earlier( k, nk, ni, np_to_p, false );
    bool was_n_blocked_by_k = point_earlier(nk,  k, ni, np_to_p, true  );
    std::swap(g.p,old_p);
    std::swap(g.time,old_time);
    std::swap(g.is_covered,old_covered);
    bool is_k_blocked_by_n  = !g.is_covered && point_earlier( k, nk, ni, np_to_p, false );
    bool is_n_blocked_by_k  = !g.is_covered && point_earlier(nk,  k, ni, np_to_p, true  );

    // debug
    if (do_extra_checks)
    {
      if (was_n_blocked_by_k)
      {
        assert( ng.my_earlier_squares.find(k) != ng.my_earlier_squares.end() );
      }
      if (was_k_blocked_by_n)
      {
        assert( g.my_earlier_squares.find(nk) != g.my_earlier_squares.end() );
      }
    }

    if (was_n_blocked_by_k && !is_n_blocked_by_k)
    {
      if (do_extra_checks)
      {
        assert( ng.my_earlier_squares.find(k) != ng.my_earlier_squares.end() );
        ng.my_earlier_squares.erase(k);
      }

      if (--ng.num_earlier == 0)
      {
        assert( ng.num_earlier >= 0 );
        potential_early_squares.push_back(nk);
      }
    }
    else if (!was_n_blocked_by_k && is_n_blocked_by_k)
    {
      if (do_extra_checks)
      {
        assert( ng.my_earlier_squares.find(k) == ng.my_earlier_squares.end() );
        ng.my_earlier_squares.insert(k);
      }

      ng.num_earlier++;
    }

    if (!g.is_covered)
    {
      if (was_k_blocked_by_n && !is_k_blocked_by_n)
      {
        if (do_extra_checks)
        {
          assert( g.my_earlier_squares.find(nk) != g.my_earlier_squares.end() );
          g.my_earlier_squares.erase(nk);
        }

        g.num_earlier--;
        assert( g.num_earlier >= 0 );
      }
      else if (!was_k_blocked_by_n && is_k_blocked_by_n)
      {
        if (do_extra_checks)
        {
          assert( g.my_earlier_squares.find(nk) == g.my_earlier_squares.end() );
          g.my_earlier_squares.insert(nk);
        }
        g.num_earlier++;
      }
    }
  }
  if (g.num_earlier==0 && !g.is_covered)
  {
    potential_early_squares.push_back(k);
  }
}

bool Grid::point_earlier( int k, int nk, int ni, std::vector<Point> &np_to_p, bool flip_sign ) const
{
  auto &g = grid_squares[k];
  auto &n = grid_squares[nk];
  // if g is earlier by time than n, then n doesn't block it
  // if n is already covered, then it doesn't block
  if (g.time<n.time || (g.time==n.time && k<nk) || n.is_covered)
  {
    return false;
  }
  //  g is later than n, and n has a candidate sample
  
  // the disk always intersects the first 8 neighbors, so is always blocked by them
  if (periodic && ni < 8)
  {
    return true;
  }
  
  // if g's sample disk doesn't intersect n's square, then n doesn't block g despite being earlier
  int nki, nkj;
  OneDtoTwo(nk, nki, nkj);
  
  if (np_to_p.empty())
  {
    return disk_intersects_square( g.p, nki, nkj );
  }

  Point p = g.p;
  if (flip_sign)
    p += np_to_p[ni];
  else
    p -= np_to_p[ni];
  return disk_intersects_square( p, nki, nkj );
}

void Grid::find_point_early(std::vector<int> &early_squares)
{
  // aggressive point-based neighbors finds more early squares
  // but checking them takes longer
  std::vector<Point> np_to_p;
  std::vector<int> ks;
  for (int gk=0; gk<(int) grid_squares.size(); ++gk)
  {
    GridSquare &g = grid_squares[gk];
    if (g.is_covered)
    {
      continue;
    }
    neighbors(gk, ks, &np_to_p);
    bool early = true;
    for (int ni=0; ni < (int) ks.size(); ++ni)
    {
      int nk = ks[ni];
      if (point_earlier(gk, nk, ni, np_to_p, false))
      {
        early = false;
        break;
      }
    }
    if (early)
    {
      early_squares.push_back(gk);
      g.is_early = true;
    }
  }
}

void Grid::count_point_early(std::vector<int> &early_squares)
{
  // aggressive point-based neighbors finds more early squares
  // but checking them takes longer
  std::vector<Point> np_to_p;
  std::vector<int> ks;
  for (int gk=0; gk<(int) grid_squares.size(); ++gk)
  {
    GridSquare &g = grid_squares[gk];
    if (g.is_covered)
    {
      continue;
    }
    neighbors(gk, ks, &np_to_p);
    for (int ni=0; ni < (int) ks.size(); ++ni)
    {
      int nk = ks[ni];
      if (point_earlier(gk, nk, ni, np_to_p, false))
      {
        g.num_earlier++;
        if (do_extra_checks)
        {
          assert( g.my_earlier_squares.find(nk) == g.my_earlier_squares.end() );
          g.my_earlier_squares.insert(nk);
        }
      }
    }
  }
  for (int gk=0; gk<(int) grid_squares.size(); ++gk)
  {
    GridSquare &g = grid_squares[gk];
    if (g.num_earlier==0 && !g.is_covered)
    {
      assert( !g.is_early );
      g.is_early = true;
      early_squares.push_back(gk);
    }
  }
}


void Grid::find_early(std::vector<int> &early_squares)
{
  std::vector<int> ks;
  for (int gk=0; gk<(int) grid_squares.size(); ++gk)
  {
    GridSquare &g = grid_squares[gk];
    if (g.is_covered || g.num_earlier>0)
    {
      continue;
    }
    neighbors(gk, ks);
    for (auto nk : ks)
    {
      auto &n_square = grid_squares[nk];
      if (n_square.is_covered)
      {
        continue;
      }
      const bool g_is_first =  g.time < n_square.time || (g.time == n_square.time && nk < gk);
      if (g_is_first)
      {
        n_square.num_earlier++;
      }
      else
      {
        g.num_earlier++;
        break;
      }
    }
    if (g.num_earlier==0)
    {
      early_squares.push_back(gk);
      g.is_early = true;
    }
  }
}

void Grid::zero_early()
{
  for (auto &g : grid_squares)
  {
    g.num_earlier=0;
  }
}

void Grid::count_early(std::vector<int> &early_squares)
{
  // count num_earlier
  {
    std::vector<int> ks;
    for (int gk=0; gk<(int) grid_squares.size(); ++gk)
    {
      GridSquare &g = grid_squares[gk];
      if (g.is_covered)
      {
        continue;
      }
      neighbors(gk, ks);
      for (auto nk : ks)
      {
        // avoid counting twice
        if (nk<=gk)
        {
          continue;
        }
        auto &n_square = grid_squares[nk];
        if (n_square.is_covered)
        {
          continue;
        }
        const bool g_is_first = g.time <= n_square.time;
        if (g_is_first)
        {
          n_square.num_earlier++;
          if (do_extra_checks)
          {
            assert( n_square.my_earlier_squares.find(gk) == n_square.my_earlier_squares.end() );
            n_square.my_earlier_squares.insert(gk);
          }
        }
        else
        {
          g.num_earlier++;
          if (do_extra_checks)
          {
            assert( g.my_earlier_squares.find(nk) == g.my_earlier_squares.end() );
            g.my_earlier_squares.insert(nk);
          }
        }
      }
    }
  }
  // find the ones with no earlier neighbors
  for (int gk=0; gk<(int) grid_squares.size(); ++gk)
  {
    GridSquare &g = grid_squares[gk];
    if (g.num_earlier==0 && !g.is_covered)
    {
      assert(!g.is_early);
      early_squares.push_back(gk);
      g.is_early = true;
    }
  }
}

// there shouldn't be any not-early or early points left, all should be accepted or empty
// this is a double-check on the number of samples we explicitly accepted
bool grid_finished(const Grid&grid, long &num_samples)
{
  num_samples = 0;
  bool ret_val = true;
  for (int i=0; i<grid.m(); ++i)
  {
    for (int j=0; j<grid.n(); ++j)
    {
      auto &g = grid.const_grid_square(i,j);
      if (g.accepted())
      {
        ++num_samples;
        continue;
      }
      if (g.empty())
        continue;;
      if (g.early() || g.not_early())
      {
        std::cout << "square (" << i << "," << j << ") is not finished" << std::endl;
        ret_val = false;
      }
    }
  }
  return ret_val;
}

void square_containing_point(Point &p, int m, int n, bool periodic, int&i , int&j)
{
  const auto a = sqrt(0.5);
  i = floor(p.x / a);
  j = floor(p.y / a);
  if (periodic)
  {
    if (i<0)
    {
      i+= m;
      p.x += m*a;
    }
    if (i>=m)
    {
      i-=m;
      p.x -= m*a;
    }
    if (j<0)
    {
      j+=n;
      p.y+=n*a;
    }
    if (j>=n)
    {
      j-=n;
      p.y-=n*a;
    }
  }
}

bool Grid::closest_other_disk(int k0, int k1, Point &p, double &dist2 ) const
{
  // find grid square of p
  int pi,pj;
  square_containing_point(p, gridm, gridn, periodic, pi,pj);
  if ((pi<0 || pi>=gridm) || (pj<0 || pj>=gridn))
  {
    // out of domain, don't care
    return false;
  }
  double pk = TwoDtoOneD(pi,pj);
  std::vector<int> ks;
  std::vector<Point> np_to_p;
  neighbors(pk, ks, &np_to_p);
  // add the intersection point's square as well
  ks.push_back(pk);
  if (!np_to_p.empty())
  {
    np_to_p.emplace_back(0.,0.);
  }
  dist2 = std::numeric_limits<double>::max();
  for (int ki = 0; ki < (int)ks.size(); ++ki)
  {
    auto nk=ks[ki];
    if (nk == k0 || nk == k1)
    {
      continue;
    }
    auto &ng = grid_squares[nk];
    if (ng.accepted())
    {
      Point v = p - ng.p;
      if (!np_to_p.empty())
        v -= np_to_p[ki];
      dist2 = std::min( dist2, normsquared(v) );
    }
  }
  if (dist2 > 1+1.e-6)
  {
    std::cout << "Disk intersection point ";
    print_point(p);
    std::cout << " in square (" << pi << "," << pj << ") is uncovered, " <<
    "too far away from a third sample, at distance " << sqrt(dist2) << std::endl;
  }
  return true;
}

bool Grid::grid_covered(bool do_progress) const
{
  // for all squares with samples
  //   for all neighbor samples
  //     check that radical axis points at distance 1 are covered by some third disk
  //     or are outside the non-periodic grid.
  // zzyk to do
  // for all boundary squares
  //   check that boundary edges are covered by disks
  //     by finding the disk closest to an endpoint, then finding the next point uncovered and recursing.

  bool ret_val = true;
  std::vector<int> ks;
  std::vector<Point> np_to_p;
  double farthest_d2 = 0.;
  for (int i=0; i<m(); ++i)
  {
    for (int j=0; j<n(); ++j)
    {
      int k = TwoDtoOneD(i,j);
      auto &g = grid_squares[k];
      if (g.accepted())
      {
        Point disk(g.p);
        neighbors(k, ks, &np_to_p);
        for (int ki = 0; ki < ks.size(); ++ki)
        {
          int nk = ks[ki];
          if (nk<k) // redundant
            continue;

          auto &ng = grid_squares[nk];
          if (ng.accepted())
          {
            Point ndisk(ng.p);
            if (!np_to_p.empty())
            {
              ndisk += np_to_p[ki];
            }
            bool disjoint;
            Point n0, n1;
            intersect_circle(disk, ndisk, disjoint, n0, n1);
            if (disjoint)
            {
              continue;
            }
            double d02, d12;
            if (closest_other_disk(k, nk, n0, d02))
              farthest_d2 = std::max( farthest_d2, d02 );
            if (closest_other_disk(k, nk, n1, d12))
              farthest_d2 = std::max( farthest_d2, d02 );
          }
        }
      }
    }
  }
  if (do_progress)
  {
    std::cout << "Farthest disk-disk intersection point to the next closest third disk is " << sqrt( farthest_d2 )
    << "; it should be <= 1 within roundoff." << std::endl;
  }
  return ret_val;
}

bool Grid::samples_separated(bool do_progress) const
{
  // for all squares
  //  if has_sample
  //  check neighbor sample distances

  bool ret_val = true;
  std::vector<int> ks;
  std::vector<Point> np_to_p;
  double closest_d2 = std::numeric_limits<double>::max();
  double maxdist = 0.; // verify the translation is working correctly.
  for (int i=0; i<m(); ++i)
  {
    for (int j=0; j<n(); ++j)
    {
      int k = TwoDtoOneD(i,j);
      auto &g = grid_squares[k];
      if (g.accepted())
      {
        neighbors(k, ks, &np_to_p);
        for (int ki = 0; ki < ks.size(); ++ki)
        {
          int nk = ks[ki];
          if (nk<k) // redundant
            continue;

          auto &ng = grid_squares[nk];
          if (ng.accepted())
          {
            Point v = g.p - ng.p;
            if (!np_to_p.empty())
            {
              v -= np_to_p[ki];
            }
            auto d2 = normsquared(v);
            closest_d2 = std::min( closest_d2, d2 );
            maxdist = std::max( maxdist, d2 );
            if (d2 < 1. - 1.e-7)
            {
              int ni, nj;
              OneDtoTwo(nk,ni,nj);

              std::cout << "Too-close samples at grid (" << i << "," << j << ") and (" << ni << "," << nj << "): ";
              print_point(g.p);
              std::cout << " and ";
              Point np = ng.p;
              if (!np_to_p.empty())
              {
                np -= np_to_p[ki];
              }
              print_point(np);
              double d = sqrt(d2);
              std::cout << " distance " << d << std::endl;
              ret_val = false;
            }
          }
        }
      }
    }
  }
  maxdist = sqrt(maxdist);
  const double theory_maxdist = sqrt( 2 + 4.5 );
  if ( maxdist > theory_maxdist + 2 * __DBL_EPSILON__)
  {
    std::cout << "Warning: max template-neighbor distance observed = " << maxdist <<
    " which is greater than the theoretical max = " << theory_maxdist << std::endl;
  }
  if (do_progress)
  {
    std::cout << "closest sample center distance is " << sqrt( closest_d2 )
    << "; it should be >= 1 within roundoff." << std::endl;
  }
  return ret_val;
}

long Grid::prepass_early(std::vector<int> &early_squares)
{
  auto startcputime = std::clock();
  long start_num_samples = (long) early_squares.size();
  find_early(early_squares);
  for (long ki = start_num_samples; ki <(long) early_squares.size(); ++ki)
  {
    accept_sample_no_counts( early_squares[ki] );
  }
  long num_accepted = (long) early_squares.size() - start_num_samples;
  if (do_iter_timing)
  {
    double cpu_duration = (std::clock() - startcputime) / (double)CLOCKS_PER_SEC;
    printf("prepass_early %2d, num_samples %8ld, cpu %5.3f, samp/sec %5.3f\n", -1, num_accepted, cpu_duration, num_accepted / cpu_duration);
    printf("%2d %8ld %5.3f %5.3f\n", -1, num_accepted, cpu_duration, num_accepted / cpu_duration);
  }
  if (do_progress)
  {
    std::cout << "prepass_early added " << num_accepted << " samples, total = " << early_squares.size() << std::endl;
  }
  return num_accepted;
}

long Grid::prepass_point_early(std::vector<int> &early_squares, int pass)
{
  auto startcputime = std::clock();
  long start_num_samples = (long) early_squares.size();
  find_point_early(early_squares);
  if (do_plot)
  {
    plot_grid("prepass_point_early", pass);
    plot_dots("prepass_point_early", pass);
  }
  for (long ki = start_num_samples; ki <(long) early_squares.size(); ++ki)
  {
    accept_sample_no_counts( early_squares[ki] );
  }
  long num_accepted = (long) early_squares.size() - start_num_samples;
  if (do_plot)
  {
    plot_grid("prepass_point_accepted", pass);
    plot_dots("prepass_point_accepted", pass);
  }
  if (do_iter_timing)
  {
    double cpu_duration = (std::clock() - startcputime) / (double)CLOCKS_PER_SEC;
    printf("\nprepass_point_early %2d, num_samples %8ld, cpu %5.3f, samp/sec %5.3f\n", -1, num_accepted, cpu_duration, num_accepted / cpu_duration);
    printf("%2d %8ld %5.3f %5.3f\n", -1, num_accepted, cpu_duration, num_accepted / cpu_duration);
  }
  if (do_progress)
  {
    std::cout << "prepass_point_early added " << num_accepted << " samples, total = " << early_squares.size() << std::endl;
  }
  return num_accepted;
}


//void Grid::report_resample_stats()
//{
//  std::vector<long>histogram(2,0);
//  for (auto &g : grid_squares)
//  {
//    if (g.times_resampled>=histogram.size())
//    {
//      histogram.resize(g.times_resampled+1, 0);
//    }
//    histogram[ g.times_resampled ]++;
//    g.times_resampled=0;
//  }
//  std::cout << "#resamples  #cells\n";
//  for (size_t i = 0; i < histogram.size(); ++i)
//  {
//    std::cout << i << "  " << histogram[i] << "\n";
//  }
//}


long Grid::sweep_grid()
{
  std::clock_t sweepstartcputime = std::clock();

  std::vector<int> early_squares;
  early_squares.reserve(grid_squares.size());
  
  AcceptWorkspace accept_workspace;
  
  long num_accepted = 0;

  // experiment with pre-passes that speed things up

  // Each prepass is O(n) time, because all squares are visited to recalculae num_earlier, even if covered already.
  
  // if the prepasses finish the grid, then we can skip the incremental part
  bool do_incremental = true;
  
  // point-neighbor prepasses. 7 is good.
  for (int pp = 0; pp<point_prepasses; ++pp)
  {
    if (do_progress)
    {
      std::cout << "\npoint-neighbor prepass " << pp << " of " << point_prepasses-1 << "\n";
    }
    long num_new = prepass_point_early(early_squares, pp);
    num_accepted += num_new;
    // report_resample_stats();
    if (do_extra_checks)
    {
      check_status();
    }
    if (num_new==0)
    {
      do_incremental = false;
      break;
    }
  }
  
  // square-neighbor prepasses. 0, because it's not as good as point-neighbors
  for (int pp = 0; pp<square_prepasses; ++pp)
  {
    if (do_progress)
    {
      std::cout << "\nsquare-neighbor prepass " << pp << " of " << square_prepasses-1 << "\n";
    }
    if (do_plot)
    {
      plot_grid("sweep_prepass_", pp);
      plot_dots("sweep_prepass_", 0);
    }

    long num_new = prepass_early(early_squares);
    num_accepted += num_new;
    zero_early();
    if (do_extra_checks)
    {
      check_status();
    }
    if (num_new==0)
    {
      do_incremental = false;
      break;
    }
  }

  if (do_incremental)
  {
    // start of main algorithm of paper
    if (do_progress)
    {
      std::cout << "\ninitial main, neighbor type = " << (do_point_neighbors ? "point" : "square" ) << "\n";
    }
    
    std::clock_t startcputime = std::clock();
    long start_num_samples = (long) early_squares.size();
    int ki = (int) start_num_samples;
    
    // find initially-early squares
    if (do_point_neighbors)
      count_point_early(early_squares);
    else
      count_early(early_squares);
    
    // debug
    if (do_extra_checks)
    {
      check_status();
    }
    if (do_print)
    {
      print_times();
      print_early_count();
      print_early();
    }
    if (do_plot)
    {
      plot_grid("sweep_", 0);
      plot_dots("sweep_", 0);
    }
    if (do_early_stats)
    {
      write_early(0);
    }
    
    // accept samples in early squares,
    //   resample neighbors
    //   update early squares based on new times for resampled squares
    int iter=0;
    long next_iter_k = ((long) early_squares.size()) -1;
    for (; ki < early_squares.size(); ++ki)
    {
      auto k = early_squares[ki];
      
      if (do_extra_checks)
      {
        auto &g = grid_squares[k];
        assert( g.num_earlier == 0 );
        assert( !g.empty() );
        assert( !g.has_sample );
        // expensive double-check for verification of correctness
        assert( do_point_neighbors ? locally_point_early(k) : locally_early(k));
      }
      
      // accept sample
      if (do_point_neighbors)
        accept_sample_update_pointneighbors(k, early_squares);
      else
        accept_sample(k, early_squares); // _update_squareneighbors
      ++num_accepted;
      
      // periodically report progress, debug
      if (ki==next_iter_k)
      {
        ++iter;
        if (do_iter_timing)
        {
          double cpu_duration = (std::clock() - startcputime) / (double)CLOCKS_PER_SEC;
          // printf("%2d %8ld %8ld   %10.3e %7.0f %6.0f  %5.3f\n", t, ag[t], num_samples, cpu_duration, ag[t] / cpu_duration, num_samples / cpu_duration, (double)ag[t]/ (double)num_samples );
          long num_iter_samples = ki - start_num_samples;
          printf("iter %2d, num_samples %8ld, cpu %5.3f, samp/sec %5.3f\n", iter, num_iter_samples, cpu_duration, num_iter_samples / cpu_duration);
          printf("%2d %8ld %5.3f %5.3f\n", iter, num_iter_samples, cpu_duration, num_iter_samples / cpu_duration);
        }
        if (do_extra_checks)
        {
          check_status();
        }
        next_iter_k = ((int) early_squares.size()) -1;
        if (do_progress)
        {
          std::cout << "iter " << iter << ", samples = " << ki << std::endl;
        }
        if (do_print)
        {
          print_times();
          print_early_count();
          print_early();
        }
        if (do_plot)
        {
          plot_grid("sweep_", iter);
          plot_dots("sweep_", iter);
        }
        if (do_early_stats)
        {
          write_early(iter);
        }
        if (do_iter_timing)
        {
          startcputime = std::clock();
          start_num_samples = ki;
        }
      }
    }
  }
    
  // i/o
  if (do_iter_timing)
  {
    double cpu_duration = (std::clock() - sweepstartcputime) / (double)CLOCKS_PER_SEC;
    printf("overall, num_samples %8ld, cpu %5.3f, samp/sec %5.3f\n", num_accepted, cpu_duration, num_accepted / cpu_duration);
    printf("   %8ld %5.3f %5.3f\n", num_accepted, cpu_duration, num_accepted / cpu_duration);
  }
  if (do_counts || do_progress)
  {
    std::cout << "number of accepted samples = " << num_accepted << "\n";
  }
  if (do_PSA_points)
  {
    print_PSA_points("points", num_accepted);
  }

  // verification of termination
  if (do_extra_checks)
  {
    long num_samples=0;
    const bool finished = grid_finished(*this, num_samples);
    if (do_progress)
    {
      std::cout << "number of grid samples = " << num_samples << std::endl;
    }
    if (do_extra_checks)
    {
      assert(num_accepted == num_samples);
      assert(finished);
    }
  }
  
  // verification of output distance properties
  if (do_extra_checks)
  {
    bool separated = samples_separated(do_progress);
    bool covered = grid_covered(do_progress);
    assert(separated);
    assert(covered);
  }

  return num_accepted;
}

// ============================================== i/o ==========

double compute_scale(int gridm, int gridn, double &dotsize)
{
  const double scale = sqrt(2) * std::min( 570. / (2+gridm), 750. / (2+gridn));
  dotsize = 0.05 * (1. - 6.0 / std::max( {8, gridm, gridn} ));
  return scale;
}


void plot_squares(int gridm, int gridn, std::ofstream &fout)
{
  const double a = sqrt(0.5);
  
  // vertical
  fout << "\n% vertical grid lines " << gridm+1 << "\n";
  fout << "newpath\n";
  for (double x = 0.; x <= gridm; ++x)
  {
    fout << x*a << " 0  moveto  ";
    fout << x*a << " " << gridn*a << " lineto\n";
  }
  fout << "stroke\n\n";

  // horizontal
  fout << "\n% horizontal grid lines " << gridn+1 << "\n";
  fout << "newpath\n";
  for (double y = 0.; y <= gridn; ++y)
  {
    fout << " 0 " << y*a << " moveto  ";
    fout << " " << gridm*a << " " << y*a << " lineto\n";
  }
  fout << "stroke\n\n";

}

void plot_grid_boundary(int gridm, int gridn, std::ofstream &fout)
{
  const double a = sqrt(0.5);
  
  // vertical
  fout << "\n% grid boundary\n";
  fout << "newpath\n";
  fout << "  0 0 moveto\n";
  fout << "  0 " << gridn*a << " lineto\n";
  fout << "  " << gridm*a << " " << gridn*a << " lineto\n";
  fout << "  " << gridm*a << " " << 0 << " lineto\n";
  fout << "closepath stroke\n";
}

void Grid::check_status() const
{
  for (int i  = 0; i < gridm; ++i)
  {
    for (int j  = 0; j < gridn; ++j)
    {
      const auto &g = const_grid_square(i,j);
      const bool accepted = g.accepted();
      const bool early = g.early();
      const bool not_early = g.not_early();
      const bool empty = g.empty();
      
      // exactly 1 should be true
      int n_true = ((int) accepted) + ((int) early) + ((int) not_early) + ((int) empty);
      if (n_true != 1)
      {
        std::cout << "error, point statuses don't add up for (" << i << "," << j << ")\n";
      }
      assert( n_true == 1);
    }
  }
}


void fout_samp(int i, int j, const Point &p, double time, bool do_circle, bool do_plot_sampletime, double dotsize, double dotfactor, std::string preamble, std::ofstream &fout)
{
  fout << "\n% " << preamble << " sample point (" << i << "," << j << ") at time " << time << "\n";
  fout << "newpath " << p.x << " " << p.y << " " << dotsize*dotfactor << " 0 360 arc closepath fill\n";
  if (do_circle)
  {
    fout << "newpath " << p.x << " " << p.y << " 1.0 0 360 arc closepath stroke\n";
  }
  if (do_plot_sampletime)
  {
    fout << "\n% sample time 'c'\n";
    fout << "/Times-Bold findfont " << 10*dotsize << " scalefont setfont\n";
    const auto a = sqrt(0.5);
    fout << (i+0.5)*a << " " << (j+0.5)*a + dotsize << " moveto (" << std::setprecision(4) << time << ") show\n";
  }
}

bool ghost_x( int m, int /*n*/, const Point &p, Point &ghost_p )
{
  const double a = sqrt(0.5);
  const double maxx = m*a;
  // const double maxy = n*a;
  const bool closex = (p.x <= 1. || p.x >= maxx-1. );
  // const bool closey = (p.y <= 1. || p.y >= maxy-1. );
  
  if (closex)
  {
    ghost_p = p;
    if (ghost_p.x <= 1.)
      ghost_p.x += maxx;
    else if (ghost_p.x >= maxx-1.)
      ghost_p.x -= maxx;
    return  true;
  }
  return false;
}

bool ghost_y( int /*m*/, int n, const Point &p, Point &ghost_p )
{
  const double a = sqrt(0.5);
  // const double maxx = m*a;
  const double maxy = n*a;
  // const bool closex = (p.x <= 1. || p.x >= maxx-1. );
  const bool closey = (p.y <= 1. || p.y >= maxy-1. );
  
  if (closey)
  {
    ghost_p = p;
    if (ghost_p.y <= 1.)
      ghost_p.y += maxy;
    else if (ghost_p.y >= maxy-1.)
      ghost_p.y -= maxy;
    return  true;
  }
  return false;
}

bool ghost_xy( int m, int n, const Point &p, Point &ghost_p )
{
  const double a = sqrt(0.5);
  const double maxx = m*a;
  const double maxy = n*a;
  const bool closex = (p.x <= 1. || p.x >= maxx-1. );
  const bool closey = (p.y <= 1. || p.y >= maxy-1. );
  
  if (closex && closey)
  {
    ghost_p = p;
    if (ghost_p.x <= 1.)
      ghost_p.x += maxx;
    else if (ghost_p.x >= maxx-1.)
      ghost_p.x -= maxx;
    if (ghost_p.y <= 1.)
      ghost_p.y += maxy;
    else if (ghost_p.y >= maxy-1.)
      ghost_p.y -= maxy;
    return  true;
  }
  return false;
}

int Grid::plot_samples_in_color(StatusFn status_fn, std::string rgb_color, std::string rgb_color_periodic,
                                bool do_circle, std::ofstream &fout, double dotsize, double linefactor) const
{
  const double default_linewidth = 0.0005;
  double lf = linefactor*plot_linewidth_factor;
  if (lf!=1.)
  {
      fout << default_linewidth*lf  << " setlinewidth\n";
  }

  int count = 0;
  fout << rgb_color << " setrgbcolor\n";
  for (int i  = 0; i < gridm; ++i)
  {
    for (int j  = 0; j < gridn; ++j)
    {
      const auto &g = const_grid_square(i,j);
      if ((g.*status_fn)())
      {
        ++count;
        fout_samp(i, j, g.p, g.time, do_circle, do_plot_sampletime, dotsize, plot_dotsize_factor, " ", fout);
      }
    }
  }
  
  if (periodic)
  {
    const double periodic_dotsize_factor = (plot_dotsize_factor >= 1.) ? (1.+plot_dotsize_factor)/2 : plot_dotsize_factor;
    if (lf!=1.)
    {
      const double periodic_linewidth_factor = (lf >= 1.) ? (1.+lf)/2 : lf;
      fout << default_linewidth*periodic_linewidth_factor  << " setlinewidth\n";
    }

    fout << rgb_color_periodic << " setrgbcolor\n";

    Point ghost_p;
    for (int i  = 0; i < gridm; ++i)
    {
      for (int j  = 0; j < gridn; ++j)
      {
        const auto &g = const_grid_square(i,j);
        if ((g.*status_fn)())
        {
          if ( ghost_x(m(),n(), g.p, ghost_p ) )
          {
            fout_samp(i, j, ghost_p, g.time, do_circle, do_plot_sampletime, dotsize, periodic_dotsize_factor, "periodic-x ghost ", fout);
          }
          if ( ghost_y(m(),n(), g.p, ghost_p ) )
          {
            fout_samp(i, j, ghost_p, g.time, do_circle, do_plot_sampletime, dotsize, periodic_dotsize_factor, "periodic-x ghost ", fout);
          }
          if ( ghost_xy(m(),n(), g.p, ghost_p ) )
          {
            fout_samp(i, j, ghost_p, g.time, do_circle, do_plot_sampletime, dotsize, periodic_dotsize_factor, "periodic-xy ghost ", fout);
          }
        }
      }
    }
  }
  if (lf!=1.)
  {
    fout << default_linewidth << " setlinewidth\n";
  }

  return count;
}

void Grid::plot_samples(std::ofstream &fout, double dotsize) const
{
  fout << "\n% samples: black=accepted, green=early, red=not early, blank = empty-polygon square \n";

  fout << "\n% not-early samples\n";
  std::string red("1.00 0.00 0.00");
  std::string lightred("2.00 0.00 0.00");
  int n_waiting = plot_samples_in_color( &GridSquare::not_early, red, lightred, false, fout, dotsize);

  fout << "\n% early samples\n";
//  std::string green("0.0 0.750 0.0");
//  std::string lightgreen("0.0 1.250 0.0");
  std::string green("0.0 0.9 0.20");
  std::string lightgreen("0.0 1.8 0.4");
  int n_early =  plot_samples_in_color( &GridSquare::early, green, lightgreen, true, fout, dotsize, 1.5);

  fout << "\n% accepted samples\n";
  std::string black("0.00 0.00 0.00");
  std::string grey("0.50 0.50 0.50");
  int n_accepted = plot_samples_in_color( &GridSquare::accepted, black, grey, true, fout, dotsize);
  
  fout << "\n% " << n_accepted << " accepted, " << n_early << " early, " << n_waiting << " not-early\n";

  // restore black color
  fout << black << "\n";
}

void grid_preamble(double scale, std::ofstream &fout)
{
  fout.precision(std::numeric_limits<double>::max_digits10);
  fout << "%!PS\n"
  << "0.000 0.000 0.000 setrgbcolor\n" // black
  << scale << " " << scale << " scale\n"
  << "0.707 0.707 translate\n"
  << "0.0005 setlinewidth\n";
}

void Grid::plot_grid(std::string fname, int id) const
{
  fname += std::to_string(id) + "-grid.ps";
  std::cout << "plotting " << fname  << std::endl;
  std::ofstream fout;
  fout.open(fname);
    
  // page scale
  double dotsize;
  const double scale = compute_scale(gridm, gridn, dotsize);

  grid_preamble(scale, fout);

  plot_squares(gridm, gridn, fout);
  plot_samples(fout, dotsize);

  fout << "\nshowpage\n";
  fout.close();
}

void Grid::plot_dots(std::string fname, int id) const
{
  fname += std::to_string(id) + "-dots.ps";
  std::cout << "plotting " << fname  << std::endl;
  std::ofstream fout;
  fout.open(fname);
    
  double dotsize;
  const double scale = compute_scale(gridm, gridn, dotsize);

  grid_preamble(scale, fout);

  plot_grid_boundary(gridm, gridn, fout);

  fout << "\n% accepted samples\n";
  std::string black("0.00 0.00 0.00");
  std::string grey("0.50 0.50 0.50");
  int ndots = plot_samples_in_color( &GridSquare::accepted, black, grey, false, fout, dotsize);
  fout << "\n% " << ndots << " accepted samples\n";

  fout << "\nshowpage\n";
  fout.close();

}

void plot_g(const Point &p, double time, const std::vector<Polygon> &polys, double scale, double dotsize, int i, int j, std::string fname)
{
  fname += "-gridsquare-polys" + std::to_string(i) + "_" + std::to_string(j) + ".ps";
  std::cout << "plotting " << fname  << std::endl;
  std::ofstream fout;
  fout.open(fname);
    
  // scale
  grid_preamble(scale, fout);

  // square
  const double a=sqrt(0.5);
  const double x0 = i*a;
  const double y0 = j*a;
  const double x1 = x0+a;
  const double y1 = y0+a;
  
  fout << "\n% grid square, unexpanded (" << i << "," << j << ") \n";
  fout << "newpath\n";
  fout << x0 << " " << y0 << " moveto\n";
  fout << x0 << " " << y1 << " lineto\n";
  fout << x1 << " " << y1 << " lineto\n";
  fout << x1 << " " << y0 << " lineto\n";
  fout << "closepath stroke\n";
  
  // polys
  for (auto &poly : polys)
  {
    // draw each segment of the polygon
    fout << "\n% poly with " << poly.size() << " segments\n";
    for (auto &s : poly)
    {
      s.plot_segment(fout);
    }
  }

  // sample
  fout << "\n% sample point (" << i << "," << j << ") at time " << time << "\n";
  fout << "newpath " << p.x << " " << p.y << " " << dotsize << " 0 360 arc closepath fill\n";

  fout << "\nshowpage\n";
  fout.close();
}

void Grid::plot_gridsquare(int i, int j, std::string fname) const
{
  double dotsize;
  double scale = compute_scale(gridm, gridn, dotsize);
  const GridSquare &g = const_grid_square(i,j);
  plot_g(g.p, g.time, g.polys, scale, dotsize, i, j, fname);
}

void plot_g(const Point &p, double time, double scale, double dotsize,
            int i, int j, const Chocks &chocks, Loop &loop, const Triangles &triangles, std::string fname)
{
  fname += "-gridsquare-triangles" + std::to_string(i) + "_" + std::to_string(j) + ".ps";
  std::cout << "plotting " << fname  << std::endl;
  std::ofstream fout;
  fout.open(fname);
    
  // scale
  grid_preamble(scale, fout);

  plot_chocks(chocks, fout);
  plot_triangles(loop, triangles, fout);
  
  if (chocks.empty() && triangles.empty() && !loop.empty())
  {
    plot_loop(loop, fout, true );
  }

  // sample
  fout << "\n% sample point (" << i << "," << j << ") at time " << time << "\n";
  fout << "newpath " << p.x << " " << p.y << " " << dotsize << " 0 360 arc closepath fill\n";

  fout << "\nshowpage\n";
  fout.close();
}

void Grid::plot_gridsquare(int i, int j, const Chocks &chocks, Loop &loop, const Triangles &triangles, std::string fname) const
{
  double dotsize;
  const double scale = compute_scale(gridm, gridn, dotsize);
  const GridSquare &g = const_grid_square(i,j);
  plot_g(g.p, g.time, scale, dotsize, i,j, chocks, loop, triangles, fname );
}

void Grid::print_PSA_points(std::string fname, long num_accepted)
{
  // count number of accepted samples if it wasn't passed in
  {
    num_accepted = 0;
    for (int i  = 0; i < gridm; ++i)
    {
      for (int j  = 0; j < gridn; ++j)
      {
        const auto &g = const_grid_square(i,j);
        if (g.accepted())
        {
          ++num_accepted;
        }
      }
    }
  }

  // when averaging files, PSA will truncate the input to the minimum number of points over all files.
  // So we randomly shuffle them to avoid a horizontal white line artifact
  std::vector<Point> shuffled_points;
  shuffled_points.reserve( num_accepted );
  for (int i  = 0; i < gridm; ++i)
  {
    for (int j  = 0; j < gridn; ++j)
    {
      const auto &g = const_grid_square(i,j);
      if (g.accepted())
      {
        shuffled_points.push_back(g.p);
      }
    }
  }
  long verified_num_accepted = shuffled_points.size();
  if (num_accepted<=0)
  {
    num_accepted = verified_num_accepted;
  }
  
  if (num_accepted != verified_num_accepted)
  {
    std::cout << "ERROR: in Grid::print_PSA_points, I was told there were going to be " << num_accepted << " points,"
    " but I found " << verified_num_accepted << " instead." << std::endl;
    std::cout << "Solutions: (1) re-call with the right number of samples passed in, or (2) recall with -1 passed in\n";
  }

  std::shuffle(shuffled_points.begin(), shuffled_points.end(), gen);
  
  // open file
  fname += std::to_string(do_PSA_points) + ".txt";
  std::cout << "printing PSA points to file " << fname  << std::endl;
  std::ofstream fout;
  fout.open(fname);
  fout.precision(16);
  
  // first line is the number of samples
  fout << verified_num_accepted << "\n";
  // print coordinates of accepted samples
  // scale coordinates to [0,1], for PSA
  const double scale = sqrt(2.) / (double) gridm;
  for (auto &p : shuffled_points)
  {
    fout << p.x * scale << "    \t" << p.y * scale << "\n";
  }
  
  fout.close();
}

void Grid::write_early(int iter)
{
  std::string fname = "early_statistics_" + std::to_string(gridm) + "X" + std::to_string(gridn) + ".txt";
  if (iter==0)
    std::cout << iter << " writing early statistics points to file " << fname  << std::endl;
  else
    std::cout << iter << " appending to " << fname << std::endl;

  // open file; clear contents if iter==0, else append
  std::ofstream fout;
  std::ios_base::openmode mode = std::ios_base::out;
  if (iter>0)
  {
    mode |= std::ios_base::app;
  }
  fout.open(fname,mode);
  
  if (iter==0)
  {
    fout << "iter accepted |num_earlier==0| |num_earlier==1| ... |num_earlier==20| empty\n";
  }
  
  // count how many grid squares have each "num_earlier" value
  std::vector<unsigned> earlier_histogram(23,0);
  for (int k=0; k<(int)grid_squares.size(); ++k)
  {
    auto &g = grid_squares[k];
    if (g.accepted())
    {
      ++earlier_histogram[0];
    }
    else if (g.empty())
    {
      ++earlier_histogram[22];
    }
    else
    {
      assert( g.num_earlier >=  0 );
      assert( g.num_earlier <= 20 );
      ++earlier_histogram[ 1 + g.num_earlier ];
    }
  }
  
  // io
  fout << iter;
  for (size_t i=0; i<earlier_histogram.size(); ++i)
  {
    fout << " " << earlier_histogram[i];
  }
  fout << "\n";
  
  fout.close();

}

void Grid::print_times() const
{
  const auto gridm = m();
  const auto gridn = n();
  
  std::cout << "Times:" << std::endl;
  for (int j = gridn-1; j>=0; --j)
  {
    for (int i = 0; i<gridm; ++i)
    {
      std::cout << const_grid_square(i,j).time << "    ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

void Grid::print_early() const
{
  const auto gridm = m();
  const auto gridn = n();
  
  std::cout << "Early:" << std::endl;
  for (int j = gridn-1; j>=0; --j)
  {
    for (int i = 0; i<gridm; ++i)
    {
      auto &g = const_grid_square(i,j);
      if (g.has_sample)
        std::cout << "a";
      else
        std::cout << ( g.is_early ? "E" : ".");
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

void Grid::print_early_count() const
{
  const auto gridm = m();
  const auto gridn = n();
  
  std::cout << "Early count:" << std::endl;
  for (int j = gridn-1; j>=0; --j)
  {
    for (int i = 0; i<gridm; ++i)
    {
      auto &g = const_grid_square(i,j);
      if (g.has_sample)
        std::cout << "a.";
      else if (g.is_covered)
        std::cout << "c.";
      else
        std::cout << g.num_earlier << ".";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}
