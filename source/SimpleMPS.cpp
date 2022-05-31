//
//  SimpleMPS.cpp
//  AvoidReject
//
//  Created by Mitchell, Scott A on 12/8/21.
//  Copyright © 2021 Mitchell, Scott A. All rights reserved.
//

#include "SimpleMPS.hpp"
#include <iostream>
#include <fstream>

long neighbor_calls=0;
long num_checks=0;

void SimpleMPS::make_grid(int m, int n)
{
  gridm = m;
  gridn = n;
  squares.resize(m*n);
  refinement_level = 0;
  twoR = 1;
  pool.reserve(2*squares.size()); // B=2 in paper for 2D
  pool2.reserve(2*squares.size());
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < n; ++j)
      pool.emplace_back(i,j);
}

Cell SimpleMPS::parent_cell(Cell subcell) const
{
  Cell pc (subcell.first / twoR, subcell.second / twoR);
  return pc;
}

long SimpleMPS::sweep_grid()
{
  long num_accepted = 0;
  if (do_progress)
  {
    std::cout << "SimpleMPS sweeping grid (" << gridm << "," << gridn << ") ";
    std::cout << (periodic? "periodic" : "chopped") << "\n";
  }
  if (do_plot)
  {
    plot_pool("pool", refinement_level);
  }
  while (!pool.empty())
  {
    if (do_progress)
    {
#ifdef DO_OPS
      std::cout << refinement_level << " " << num_accepted << " ";
#else
      std::cout << refinement_level << " iter, samples = " << num_accepted << ", pool size = " << pool.size() << std::endl;
#endif
    }
    
    // actual algorithm
    num_accepted += sample_grid();
    
    // debug quit after top-level sample
    // return num_accepted;
    
    refine();

    // debug quit after top-level sample
    // return num_accepted;

    // extra checks not needed in practice, rarely gets to that level, and if it does, it is just an efficiency issue
    // if (refinement_level<level_limit)
    //  refine();
    // else
    //   repool(); // without refining, just check if square==point is covered

    if (do_plot)
    {
      plot_pool("pool", refinement_level);
      plot_grid("iter", refinement_level);
    }
  }
  if (do_progress)
  {
    std::cout << refinement_level << " iter done, finished with samples = " << num_accepted << std::endl;
    std::cout << "a/twoR = " << sqrt(0.5)/twoR << " min box side size " << std::endl;
  }
  return num_accepted;
}

void SimpleMPS::random_point( const Cell &c, Point &p )
{
  const auto a = sqrt(0.5);
  const auto ell = a / (double) twoR;
  p.x = (c.first  + myrand.u()) * ell;
  p.y = (c.second + myrand.u()) * ell;
}

long SimpleMPS::sample_grid()
{
  // workspace
  std::vector<int> ks;
  std::vector<Point> np_to_p;
  
#ifdef DO_OPS
  neighbor_calls = 0;
#endif
  
  long num_samples = 0;
  const double tries_per_candidate = 0.3; // 0.5; // A=0.3 in paper for 2D. About the same speed for A=0.3--0.6
  const long max_tries = 16 + (pool.size() * tries_per_candidate);
  int t=0;
  for (; t<max_tries && !pool.empty(); ++t)
  {
    size_t ki = (size_t) floor( myrand.u()* (double)pool.size() );
    assert( ki < pool.size() );
    auto &subcell = pool[ki];
    auto pc = parent_cell(subcell);
    long pk = TwoDtoOneD(gridm, gridn, pc.first, pc.second);

    // not updating the pool appears to be faster!
    // the cost of the swap is high in terms of memory locality, and we are unlikely to pick the same subcell again.
    
    // with pool updates
    /*
    if (squares[pk].has_sample)
    {
      // remove subcell from pool
      std::swap( subcell, pool.back() );
      pool.pop_back();
      continue;
    }
    random_point( subcell, squares[pk].p);
    if (!covered_point(pk, ks, np_to_p))
    {
      // accept sample
      squares[pk].has_sample = true;
      ++num_samples;
      
      // remove subcell from pool
      std::swap( subcell, pool.back() );
      pool.pop_back();
    }
     */
    
    // without pool updates
    if (!squares[pk].has_sample)
    {
      random_point( subcell, squares[pk].p);
      if (!covered_point(pk, ks, np_to_p))
      {
        // accept sample
        squares[pk].has_sample = true;
        ++num_samples;
      }
    }

  }
  
#ifdef DO_OPS
//  std::cout << "number of neighbor calls = " << neighbor_calls << "\n";
//  std::cout << "number of darts thrown = " << t << " out of " << max_tries << "\n";
//  std::cout << "number of successful samples = " << num_samples << "\n";
  std::cout <<  neighbor_calls << " " << t << " " << max_tries << " " << num_samples << " ";
#endif
  
  return num_samples;
}

void SimpleMPS::corners(const Cell &c, Point &p00, Point &p10, Point &p01, Point &p11) const
{
  const auto a = sqrt(0.5);
  const double ell = a / (double) twoR; // side length of a subcell
  p00.x = c.first * ell;
  p00.y = c.second * ell;
  p10.x = p00.x + ell;
  p10.y = p00.y;
  p01.x = p00.x;
  p01.y = p00.y + ell;
  p11.x = p10.x;
  p11.y = p01.y;
}

/* unused
bool SimpleMPS::covered_cell(Cell c) const
{
  Point p00, p10, p01, p11;
  corners(c, p00, p10, p01, p11);
  auto pc = parent_cell(c);
  auto pk = TwoDtoOneD(gridm, gridn, pc.first, pc.second);
  if (squares[pk].has_sample)
  {
    return true;
  }
  std::vector<int> ks;
  std::vector<Point> np_to_p;
#ifdef DO_OPS
  ++neighbor_calls;
#endif
  neighbors(gridm, gridn, periodic, pk, ks, &np_to_p);
  for (int ki=0; ki<ks.size(); ++ki)
  {
    auto nk = ks[ki];
    if (squares[nk].has_sample)
    {
#ifdef DO_OPS
      ++num_checks;
#endif
      Point p = squares[nk].p;
      if (!np_to_p.empty())
      {
        p += np_to_p[ki];
      }
      if ( normsquared(p00-p) <= 1. && normsquared(p10-p) <= 1. && normsquared(p01-p) <= 1. && normsquared(p11-p) <= 1. )
      {
        return true;
      }
    }
  }
  return false;
}
*/

bool SimpleMPS::covered_point(long k, std::vector<int> &ks, std::vector<Point> &np_to_p) const
{
  assert(k>=0 && k < squares.size());
  if (squares[k].has_sample)
  {
    return true;
  }
#ifdef DO_OPS
  ++neighbor_calls;
#endif
  neighbors(gridm, gridn, periodic, (int)k, ks, &np_to_p);
  for (int ki=0; ki<(int)ks.size(); ++ki)
  {
    auto nk = ks[ki];
    // zzyk here is a conditional that might take variable time
    if (squares[nk].has_sample)
    {
      Point p = squares[nk].p;
      if (!np_to_p.empty())
      {
        p += np_to_p[ki];
      }

      if ( normsquared(p - squares[k].p) <= 1.)
      {
        return true;
      }
    }
  }
  return false;
}

//void SimpleMPS::refine()
//{
//  ++refinement_level;
//  twoR *=2;
//  pool2.clear();
//  for (auto &c : pool)
//  {
//    const auto ci2 = c.first*2;
//    const auto cj2 = c.second*2;
//    Cell cc00( ci2,   cj2);    if (!covered_cell(cc00)) pool2.push_back(cc00);
//    Cell cc10( ci2+1, cj2);    if (!covered_cell(cc10)) pool2.push_back(cc10);
//    Cell cc01( ci2,   cj2+1);  if (!covered_cell(cc01)) pool2.push_back(cc01);
//    Cell cc11( ci2+1, cj2+1);  if (!covered_cell(cc11)) pool2.push_back(cc11);
//  }
//  std::swap(pool,pool2);
//}

void SimpleMPS::refine()
{
  Point p00, p10, p20, p01, p02, p11, p12, p21, p22;
  std::vector<int> ks;
  std::vector<Point> np_to_p;

  pool2.clear();

#ifdef DO_OPS
  size_t pool2_start_capacity = pool2.capacity();
  neighbor_calls = 0;
#endif

  for (auto &c : pool)
  {
    
    auto pc = parent_cell(c);
    auto pk = TwoDtoOneD(gridm, gridn, pc.first, pc.second);
    if (squares[pk].has_sample)
    {
      continue;
    }
    corners(c, p00, p20, p02, p22);
    p10.x = (p00.x+p20.x)*0.5;
    p10.y = p00.y;
    p01.x = p00.x;
    p01.y = (p00.y+p02.y)*0.5;
    p11.x = p10.x;
    p11.y = p01.y;
    p12.x = p10.x;
    p12.y = p22.y;
    p21.x = p22.x;
    p21.y = p01.y;

    bool covered00(false), covered01(false), covered10(false), covered11(false);
    
#ifdef DO_OPS
    ++neighbor_calls;
#endif
    neighbors(gridm, gridn, periodic, pk, ks, &np_to_p);
    for (int ki=0; ki<ks.size(); ++ki)
    {
      auto nk = ks[ki];
      if (squares[nk].has_sample)
      {
        Point p = squares[nk].p;
        if (!np_to_p.empty())
        {
          p += np_to_p[ki];
        }
        // if it doesn't covere the center, it doesn't cover any square.
        // zzyk here is a conditional that might take variable time
        if (normsquared(p11-p) <= 1.)
        {
          if (normsquared(p10-p) <= 1.)
          {
            if (!covered00 && (normsquared(p00-p) <= 1.) && (normsquared(p01-p) <= 1.))
            {
              covered00 = true;
            }
            if (!covered10 && (normsquared(p20-p) <= 1.) && (normsquared(p21-p) <= 1.))
            {
              covered10 = true;
            }
          }
          if (normsquared(p12-p) <= 1.)
          {
            if ((!covered01 && normsquared(p02-p) <= 1.) && (normsquared(p01-p) <= 1.))
            {
              covered01 = true;
            }
            if ((!covered11 && normsquared(p22-p) <= 1.) && (normsquared(p21-p) <= 1.))
            {
              covered11 = true;
            }
          }
          if (covered00 && covered01 && covered10 && covered11)
            break;
        }
      }
    }
    const auto ci2 = c.first*2;
    const auto cj2 = c.second*2;
    if (!covered00) pool2.emplace_back( ci2  , cj2  );
    if (!covered01) pool2.emplace_back( ci2  , cj2+1);
    if (!covered10) pool2.emplace_back( ci2+1, cj2  );
    if (!covered11) pool2.emplace_back( ci2+1, cj2+1);
  }

#ifdef DO_OPS
  if (pool2_start_capacity != pool2.capacity())
    std::cout << "Warning, pool2 capacity was increased\n";
//  std::cout << "number of neighbor calls = " << neighbor_calls << "\n";
//  std::cout << "number of cells = " << pool.size() << "\n";
//  std::cout << "number of uncovered child cells = " << pool2.size() << "\n";
  std::cout << neighbor_calls << " " << pool.size() << " " << pool2.size() << "\n";
#endif

  ++refinement_level;
  twoR *=2;
  std::swap(pool,pool2);
}


int SimpleMPS::plot_samples(std::ofstream &fout, double dotsize, bool do_circle) const
{
  fout << "\n% simpleMPS accepted samples\n";
  std::string black("0.00 0.00 0.00");
  std::string grey("0.50 0.50 0.50");

  int count = 0;
  fout << black << " setrgbcolor\n";
  for (int i  = 0; i < gridm; ++i)
  {
    for (int j  = 0; j < gridn; ++j)
    {
      auto k = TwoDtoOneD(gridm, gridn, i,j);
      const auto &g = squares[k];
      if (g.has_sample)
      {
        ++count;
        const Point &p = g.p;
        fout << "\n% sample point (" << i << "," << j << ")\n";
        fout << "newpath " << p.x << " " << p.y << " " << dotsize << " 0 360 arc closepath fill\n";
        if (do_circle)
        {
          fout << "newpath " << p.x << " " << p.y << " 1.0 0 360 arc closepath stroke\n";
        }
      }
    }
  }
  if (periodic)
  {
    fout << grey << " setrgbcolor\n";

    for (int i  = 0; i < gridm; ++i)
    {
      for (int j  = 0; j < gridn; ++j)
      {
        auto k = TwoDtoOneD(gridm, gridn, i,j);
        const auto &g = squares[k];
        if (g.has_sample)
        {
          const double a = sqrt(0.5);
          const double maxx = gridm*a;
          const double maxy = gridn*a;
          const bool closex = (g.p.x <= 1. || g.p.x >= maxx-1. );
          const bool closey = (g.p.y <= 1. || g.p.y >= maxy-1. );
          
          if (closex)
          {
            Point p = g.p;
            if (p.x <= 1.)
              p.x += maxx;
            else if (p.x >= maxx-1.)
              p.x -= maxx;
            
            fout << "\n% periodic-x ghost sample point (" << i << "," << j << ")\n";
            fout << "newpath " << p.x << " " << p.y << " " << dotsize << " 0 360 arc closepath fill\n";
            if (do_circle)
            {
              fout << "newpath " << p.x << " " << p.y << " 1.0 0 360 arc closepath stroke\n";
            }
          }

          if (closey)
          {
            Point p = g.p;
            if (p.y <= 1.)
              p.y += maxy;
            else if (p.y >= maxy-1.)
              p.y -= maxy;
            
            fout << "\n% periodic-y ghost sample point (" << i << "," << j << ")\n";
            fout << "newpath " << p.x << " " << p.y << " " << dotsize << " 0 360 arc closepath fill\n";
            if (do_circle)
            {
              fout << "newpath " << p.x << " " << p.y << " 1.0 0 360 arc closepath stroke\n";
            }
          }

          if (closex && closey)
          {
            Point p = g.p;
            if (p.x <= 1.)
              p.x += maxx;
            else if (p.x >= maxx-1.)
              p.x -= maxx;
            if (p.y <= 1.)
              p.y += maxy;
            else if (p.y >= maxy-1.)
              p.y -= maxy;
            
            fout << "\n% periodic-xy ghost sample point (" << i << "," << j << ")\n";
            fout << "newpath " << p.x << " " << p.y << " " << dotsize << " 0 360 arc closepath fill\n";
            if (do_circle)
            {
              fout << "newpath " << p.x << " " << p.y << " 1.0 0 360 arc closepath stroke\n";
            }
          }

        }
      }
    }
    fout << black << " setrgbcolor\n";
  }
  return count;
}



void SimpleMPS::plot_grid(std::string fname, int id) const
{
  fname += std::to_string(id) + "-SimpleMPSgrid.ps";
  std::cout << "plotting " << fname  << std::endl;
  std::ofstream fout;
  fout.open(fname);
    
  // page scale
  double dotsize;
  const double scale = compute_scale(gridm, gridn, dotsize);

  grid_preamble(scale, fout);

  plot_squares(gridm, gridn, fout);
  plot_samples(fout, dotsize, true);

  fout << "\nshowpage\n";
  fout.close();
}

void SimpleMPS::plot_pool(std::string fname, int id) const
{
  fname += std::to_string(id) + "-SimpleMPSpool.ps";
  std::cout << "plotting " << fname  << std::endl;
  std::ofstream fout;
  fout.open(fname);
    
  // page scale
  double dotsize;
  const double scale = compute_scale(gridm, gridn, dotsize);

  grid_preamble(scale, fout);

  Point p00, p10, p01, p11;
  std::string black("0.00 0.00 0.00");
  std::string red("0.90 0.25 0.00");

  // shaded cells
  fout << red << " setrgbcolor\n";
  for (auto &c : pool)
  {
    fout << "\n% cell (" << c.first << "," << c.second << ") shade\n";
    corners(c, p00, p10, p01, p11);
    fout << p00.x << " " << p00.y << " moveto\n";
    fout << p10.x << " " << p10.y << " lineto\n";
    fout << p11.x << " " << p11.y << " lineto\n";
    fout << p01.x << " " << p01.y << " lineto\n";
    fout << "closepath fill\n";
  }

  // outline of cells
  fout << black << " setrgbcolor\n";
  for (auto &c : pool)
  {
    fout << "\n% cell (" << c.first << "," << c.second << ") corners\n";
    corners(c, p00, p10, p01, p11);
    fout << p00.x << " " << p00.y << " moveto\n";
    fout << p10.x << " " << p10.y << " lineto\n";
    fout << p11.x << " " << p11.y << " lineto\n";
    fout << p01.x << " " << p01.y << " lineto\n";
    fout << "closepath stroke\n";
  }
  
  plot_samples(fout, dotsize, true);

  fout << "\nshowpage\n";
  fout.close();

}

void SimpleMPS::plot_dots(std::string fname, int id) const
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
  plot_samples(fout, dotsize, false);

  fout << "\nshowpage\n";
  fout.close();

}
