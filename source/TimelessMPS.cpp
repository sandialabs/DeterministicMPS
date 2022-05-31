//
//  TimelessMPS.cpp
//  AvoidReject
//
//  Created by Mitchell, Scott A on 1/24/22.
//  Copyright © 2022 Mitchell, Scott A. All rights reserved.
//

#include "TimelessMPS.hpp"

#include <iostream>
#include <fstream>
#include <ostream>
#include <iomanip>

#include "Grid.hpp"

void TimelessGrid::update_area(int k, double new_area)
{
  assert(new_area>=0.);
  squares[k].area = new_area;

  // propagate changes up tree
  int nodei = k;
  int nodej = k%2 == 0 ? k+1 : k-1;
  int node = k/2;
  assert( nodei<squares.size() );
  assert(node <tree[0].size());
  tree[0][node] = (nodej < squares.size() ? squares[nodei].area + squares[nodej].area : squares[nodei].area);
  assert(tree[0][node] <= 1. + 1e-12);
  int layer = 1;
  while (layer<tree.size())
  {
    nodei = node;
    nodej = node%2 == 0 ? node+1 : node-1;
    node = node/2;
    assert(nodei<tree[layer-1].size());
    assert(node <tree[layer].size());
    tree[layer][node] = (nodej<tree[layer-1].size() ? tree[layer-1][nodei] + tree[layer-1][nodej] : tree[layer-1][nodei]);
    ++layer;
  }
}


// std::vector< std::vector<double> > tree;
void TimelessGrid::initialize_tree()
{
  const int gridsize = gridm*gridn;
  size_t num_layers = ceil( log2( gridsize ) );
  tree.resize(num_layers);
  int priorsize = gridsize, layersize;
  int layer = 0;
  do
  {
    assert(layer<tree.size());
    layersize = (priorsize+1)/2;
    tree[layer++].resize(layersize);
    priorsize = layersize;
  } while (layersize>1);
  assert( layer == num_layers );
  
  // optimization: faster way to do this for all nodes at once
  for (int k=0; k<gridsize; ++k)
  {
    update_area(k,0.5);
  }
}

int TimelessGrid::random_square()
{
  const auto u = myrand.u();
  auto a = tree.back()[0]*u;
  int k = treewalk(a);
  return k;
}

int TimelessGrid::treewalk(double a)
{
  int k=0;
  for (auto it = tree.rbegin()+1; it != tree.rend(); it++)
  {
    const int i = k*2, j = 2*k+1;
    if (j >= it->size())
    {
      k = i;
      continue;
    }
    const auto &ai = (*it)[i], &aj = (*it)[j];
    assert( a <= ai + aj );
    // go ? left : right
    if ( ai>0. && (a <= ai || aj<=0.) )
    {
      k = i;
    }
    else
    {
      k = j;
      a -= ai;
    }
  }
  assert(k>=0 && k< ceil(squares.size()/2.0));
  const int i = k*2, j = k*2+1;
  const auto ai = squares[i].area;
  if ( j>= squares.size() )
  {
    k = i;
  }
  else
  {
    assert( j < squares.size() );
    const auto aj = squares[j].area;
    if ( ai>0. && (a <= ai || aj<=0.) )
    {
      k = i;
    }
    else
    {
      k = j;
      a -= ai;
    }
  }
  assert( squares[k].area > 0. );
  assert( squares[k].area >= a );
  assert( ai >= 0. );
  return k;
}

// there shouldn't be any uncovered squares
// also double-check on the number of samples we explicitly accepted
bool grid_finished(const TimelessGrid&grid, long &num_samples)
{
  num_samples = 0;
  bool ret_val = true;
  for (int i=0; i<grid.m(); ++i)
  {
    for (int j=0; j<grid.n(); ++j)
    {
      auto &g = grid.const_square(i,j);
      if (g.has_sample)
      {
        assert(g.is_covered);
        ++num_samples;
        continue;
      }
      if (g.is_covered)
      {
        continue;
      }
      {
        std::cout << "square (" << i << "," << j << ") is not finished; it doesn't have a sample, but is not covered.\n" << std::endl;
        ret_val = false;
      }
    }
  }
  return ret_val;
}

bool closest_other_disk(const TimelessGrid &grid, int k0, int k1, Point &p, double &dist2 )
{
  const auto gridm = grid.m();
  const auto gridn = grid.n();
  // find grid square of p
  int pi,pj;
  square_containing_point(p, gridm, gridn, grid.periodic, pi,pj);
  if ((pi<0 || pi>=gridm) || (pj<0 || pj>=gridn))
  {
    // out of domain, don't care
    return false;
  }
  auto pk = TwoDtoOneD(gridm, gridn, pi,pj);
  std::vector<int> ks;
  std::vector<Point> np_to_p;
  neighbors(gridm,gridn,grid.periodic,pk, ks, &np_to_p);
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
    auto &ng = grid.const_square(nk);
    if (ng.has_sample)
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
bool grid_covered(const TimelessGrid&grid, bool do_progress)
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
  auto m = grid.m();
  auto n = grid.n();
  for (int i=0; i<m; ++i)
  {
    for (int j=0; j<n; ++j)
    {
      int k = TwoDtoOneD(m,n,i,j);
      auto &g = grid.const_square(i,j);
      if (g.has_sample)
      {
        Point disk(g.p);
        neighbors(m,n,grid.periodic,k,ks,&np_to_p);
        for (int ki = 0; ki < ks.size(); ++ki)
        {
          int nk = ks[ki];
          if (nk<k) // redundant
            continue;

          auto &ng = grid.const_square(nk);
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
            if (closest_other_disk(grid, k, nk, n0, d02))
              farthest_d2 = std::max( farthest_d2, d02 );
            if (closest_other_disk(grid, k, nk, n1, d12))
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

bool samples_separated(const TimelessGrid&grid, bool do_progress)
{
  // for all squares
  //  if has_sample
  //  check neighbor sample distances

  bool ret_val = true;
  std::vector<int> ks;
  std::vector<Point> np_to_p;
  double closest_d2 = std::numeric_limits<double>::max();
  double maxdist = 0.; // verify the translation is working correctly.
  auto m = grid.m();
  auto n = grid.n();
  for (int i=0; i<m; ++i)
  {
    for (int j=0; j<n; ++j)
    {
      int k = TwoDtoOneD(m,n,i,j);
      auto &g = grid.const_square(k);
      if (g.has_sample)
      {
        neighbors(m,n,grid.periodic,k,ks,&np_to_p);
        for (int ki = 0; ki < ks.size(); ++ki)
        {
          int nk = ks[ki];
          if (nk<k) // redundant
            continue;

          auto &ng = grid.const_square(nk);
          if (ng.has_sample)
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
              OneDtoTwo(m,n,nk,ni,nj);

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

void TimelessGrid::make(int m, int n)
{
  gridm = m;
  gridn = n;
  squares.resize(m*n);
  initialize_tree();
  
  // create points, as a speedup
  for (int k=0; k<(int)squares.size(); ++k)
  {
    int i, j;
    OneDtoTwo(k,i,j);
    
    auto &g = squares[k];
    g.p = random_point_in_square(i,j,myrand);
  }
  accept_workspace.ks.reserve(20);
  accept_workspace.np_to_p.reserve(20);
  scoop_workspace.reserve(4); // max of 4 voids is possible
  sample_workspace.resize(4); // max of 4 voids is possible
}


long TimelessGrid::sweep()
{
  int num_plots = do_plot || do_progress ? 10 : 1;
  const int plot_period = squares.size() / (2.88 * num_plots);
  int plot_countdown = plot_period;

  int num_samples = 0;
  if ( do_plot )
  {
    plot_grid("timeless_", num_samples);
    plot_dots("timeless_", num_samples);
  }
  double area_left;
  do
  {
    auto k = random_square();
    accept(k);
    ++num_samples;
    area_left = tree.back()[0];

    if (--plot_countdown == 0 || area_left == 0.)
    {
      if ( do_plot )
      {
        plot_grid("timeless_", num_samples);
        plot_dots("timeless_", num_samples);
      }
      if (do_progress)
      {
        // add number of non-empty squares left
        std::cout << "accepted sample " << k << std::endl;
        std::cout << num_samples << " samples, " << area_left << " area left\n";
      }
      plot_countdown = plot_period;
    }
  } while (area_left>0.);
  
  if (do_extra_checks)
  {
    long audit_num_samples;
    bool finished = grid_finished(*this, audit_num_samples);
    if (!finished)
    {
      std::cout << "I thought I finished sampling, but their are unfinished squares\n";
    }
    if (num_samples != audit_num_samples)
    {
      std::cout << "error, I thought I accepted " << num_samples << " samples, but checking the grid I found " << audit_num_samples << "\n";
    }
    bool covered = grid_covered(*this, do_progress);
    bool separated = samples_separated(*this, do_progress);
    assert(finished && num_samples == audit_num_samples);
    assert(covered);
    assert(separated);

  }
  
  return num_samples;
}

void TimelessGrid::accept(int k)
{
  auto &g = squares[k];
  assert( !g.is_covered );
  g.set_accepted();
  update_area(k,0.);
  
  auto &ks = accept_workspace.ks;
  auto &np_to_p = accept_workspace.np_to_p;
  // neighbors(k, ks, &np_to_p);
  ::neighbors(gridm, gridn, periodic, (int)k, ks, &np_to_p);
  auto disk = g.p;
  for (int ki = 0; ki < ks.size(); ++ki)
  {
    int nk = ks[ki];
    auto &n = squares[nk];
    int ni, nj;
    OneDtoTwo(nk, ni, nj);
    if (!np_to_p.empty()) disk = (g.p - np_to_p[ki]);
    if (n.is_covered || !disk_intersects_square( disk, ni, nj ))
      continue;
    double new_area;
    if ( scoop_square(nk, disk, new_area ) )
    {
      update_area(nk, new_area );
    }
  }
}

bool TimelessGrid::scoop_square(int k, Point &disk, double &new_area)
{
  auto &g = squares[k];

  assert(!g.is_covered); // should have caught this before calling
  assert(g.area > 0.); // ditto

  // build the square poly if this is the first time
  bool check_disk_covers_polygon = true;
  if (g.polys.empty())
  {
    int i, j;
    OneDtoTwo(k, i, j);
    if (::disk_covers_square(disk, i,j))
    {
      g.is_covered = true;
      new_area = 0.;
      return true;
    }
    else
    {
      check_disk_covers_polygon = false;
      const auto a = sqrt(0.5);
      g.polys.reserve(4);
      g.polys.emplace_back( square_to_poly( (double) i*a, (double) j*a ) );
    }
  }

  // scoop
  scoop_polys(g.polys, disk, scoop_workspace, check_disk_covers_polygon);
  
  // resample if not covered,
  //   the sample may not be used, but it's fast compared to generating the chocks and triangulation we needed for the area calculation
  if (g.polys.empty())
  {
    g.is_covered = true;
    new_area = 0.;
    return true;
  }
  else
  {
    // optimization: check if g.polys is unchanged
    new_area = resample(k);
    return (new_area != g.area);
  }
}

double TimelessGrid::resample(int k)
{
  auto &g = squares[k];
  // gs.times_resampled++;
  
  // shouldn't be calling this on an empty square
  assert( !g.is_covered );
  assert( !g.polys.empty() );
  
  int i, j;
  OneDtoTwo(k, i,j);
  double polys_area;
  int error_code, which_poly;
  ::sample(i, j, g.polys, myrand, sample_workspace, do_simple_triangulation,
         which_poly, g.p, polys_area, error_code);
  assert(polys_area>=0.);
  
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
    if (!check_point_in_poly(i, j, g.p, g.polys))
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
    plot_square(i,j,"sample_error");
    plot_square(i,j,chocks,loop,triangles,"sample_error");
    
    Chocks nochocks;
    Triangles notriangles;
    Loop noloop;
    plot_square(i,j, chocks, loop, notriangles, "sample_error-chocksonly");
    plot_square(i,j, nochocks, loop, notriangles, "sample_error-looponly");
    plot_square(i,j, nochocks, loop, triangles, "sample_error-trianglesonly");
  }
  return polys_area;
}

void TimelessGrid::plot_square(int i, int j, const Chocks &chocks, Loop &loop, const Triangles &triangles, std::string fname) const
{
  double dotsize;
  const double scale = compute_scale(gridm, gridn, dotsize);
  auto &g = const_square(i,j);
  plot_g(g.p, 0., scale, dotsize, i,j, chocks, loop, triangles, fname );
}

void TimelessGrid::plot_square(int i, int j, std::string fname) const
{
  double dotsize;
  double scale = compute_scale(gridm, gridn, dotsize);
  auto &g = const_square(i,j);
  plot_g(g.p, 0., g.polys, scale, dotsize, i, j, fname);
}


void TimelessGrid::plot_samples(std::ofstream &fout, double dotsize, bool accepted_dots_only) const
{
  int n_waiting = 0;
  if (!accepted_dots_only)
  {
    fout << "\n% samples: black=accepted, red=not accepted, blank = empty-polygon square \n";
    
    fout << "\n% not-accepted samples\n";
    std::string red("1.00 0.00 0.00");
    std::string lightred("2.00 0.00 0.00");
    
    bool do_circle = false;
    fout << red << " setrgbcolor\n";
    for (int i  = 0; i < gridm; ++i)
    {
      for (int j  = 0; j < gridn; ++j)
      {
        const auto &g = const_square(i,j);
        if (!g.is_covered)
        {
          ++n_waiting;
          fout_samp(i, j, g.p, 0., do_circle, false, dotsize, plot_dotsize_factor, " ", fout);
        }
      }
    }
    if (periodic)
    {
      fout << lightred << " setrgbcolor\n";
      Point ghost_p;
      const double periodic_dotsize_factor = (plot_dotsize_factor >= 1.) ? (1.+plot_dotsize_factor)/2 : plot_dotsize_factor;
      for (int i  = 0; i < gridm; ++i)
      {
        for (int j  = 0; j < gridn; ++j)
        {
          const auto &g = const_square(i,j);
          if (!g.is_covered)
          {
            if ( ghost_x(m(),n(), g.p, ghost_p ) )
            {
              fout_samp(i, j, ghost_p, 0., do_circle, false, dotsize, periodic_dotsize_factor, "periodic-x ghost ", fout);
            }
            if ( ghost_y(m(),n(), g.p, ghost_p ) )
            {
              fout_samp(i, j, ghost_p, 0., do_circle, false, dotsize, periodic_dotsize_factor, "periodic-x ghost ", fout);
            }
            if ( ghost_xy(m(),n(), g.p, ghost_p ) )
            {
              fout_samp(i, j, ghost_p, 0., do_circle, false, dotsize, periodic_dotsize_factor, "periodic-xy ghost ", fout);
            }
          }
        }
      }
    }
  }
  
  fout << "\n% accepted samples\n";
  std::string black("0.00 0.00 0.00");
  std::string grey("0.50 0.50 0.50");
  
  int n_accepted = 0;
  bool do_circle = !accepted_dots_only;
  fout << black << " setrgbcolor\n";
  for (int i  = 0; i < gridm; ++i)
  {
    for (int j  = 0; j < gridn; ++j)
    {
      const auto &g = const_square(i,j);
      if (g.has_sample)
      {
        ++n_accepted;
        fout_samp(i, j, g.p, 0., do_circle, false, dotsize, plot_dotsize_factor, " ", fout);
      }
    }
  }
  if (periodic)
  {
    fout << grey << " setrgbcolor\n";
    Point ghost_p;
    const double periodic_dotsize_factor = (plot_dotsize_factor >= 1.) ? (1.+plot_dotsize_factor)/2 : plot_dotsize_factor;
    for (int i  = 0; i < gridm; ++i)
    {
      for (int j  = 0; j < gridn; ++j)
      {
        const auto &g = const_square(i,j);
        if (g.has_sample)
        {
          if ( ghost_x(m(),n(), g.p, ghost_p ) )
          {
            fout_samp(i, j, ghost_p, 0., do_circle, false, dotsize, periodic_dotsize_factor, "periodic-x ghost ", fout);
          }
          if ( ghost_y(m(),n(), g.p, ghost_p ) )
          {
            fout_samp(i, j, ghost_p, 0., do_circle, false, dotsize, periodic_dotsize_factor, "periodic-x ghost ", fout);
          }
          if ( ghost_xy(m(),n(), g.p, ghost_p ) )
          {
            fout_samp(i, j, ghost_p, 0., do_circle, false, dotsize, periodic_dotsize_factor, "periodic-xy ghost ", fout);
          }
        }
      }
    }
  }
  
  fout << "\n% " << n_accepted << " accepted, " << n_waiting << " not-accepted placeholders\n";

  // restore black color
  fout << black << "\n";
}

void TimelessGrid::plot_grid(std::string fname, int id) const
{
  fname += std::to_string(id) + "-timelessgrid.ps";
  std::cout << "plotting " << fname  << std::endl;
  std::ofstream fout;
  fout.open(fname);
    
  // page scale
  double dotsize;
  const double scale = compute_scale(gridm, gridn, dotsize);

  grid_preamble(scale, fout);

  plot_squares(gridm, gridn, fout);
  plot_samples(fout, dotsize, false);

  fout << "\nshowpage\n";
  fout.close();
}

void TimelessGrid::plot_dots(std::string fname, int id) const
{
  fname += std::to_string(id) + "-timelessdots.ps";
  std::cout << "plotting " << fname  << std::endl;
  std::ofstream fout;
  fout.open(fname);
    
  double dotsize;
  const double scale = compute_scale(gridm, gridn, dotsize);

  grid_preamble(scale, fout);

  plot_grid_boundary(gridm, gridn, fout);
  plot_samples(fout, dotsize, true);

  fout << "\nshowpage\n";
  fout.close();

}
