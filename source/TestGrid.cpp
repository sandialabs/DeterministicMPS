//
//  TestGrid.cpp
//  AvoidReject
//
//  Created by Mitchell, Scott A on 11/12/21.
//  Copyright © 2021 Mitchell, Scott A. All rights reserved.
//

#include "TestGrid.hpp"

#include <iostream>
#include <utility>
#include <set>
#include <algorithm>
#include <ctime>
#include <chrono>

bool TestGrid::do_simple_triangulation = true;
bool TestGrid::do_early_stats = false;

void TestGrid::print_all_neighbors(const Grid &grid)
{
  int gridm = grid.m();
  int gridn = grid.n();
  
  std::vector<int> ks, is, js;
  for (int i = 0; i < gridm; ++i)
  {
    for (int j = 0; j < gridn; ++j)
    {
      // grid.neighbors(i,j, is,js);
      grid.neighbors( grid.TwoDtoOneD(i,j), ks);
      is.clear();
      js.clear();
      for (auto nk : ks)
      {
        int ni,nj;
        grid.OneDtoTwo(nk, ni,nj);
        is.push_back(ni);
        js.push_back(nj);
      }
      assert( is.size() == js.size());
      
      std::cout << "(" << i << "," << j << ") has " << is.size() << " neighbors:\n    ";
      for (int k = 0; k < is.size(); ++k)
      {
        std::cout << "(" << is[k] << "," << js[k] << ") ";
      }
      std::cout << std::endl;
    }
  }
}

void TestGrid::check_all_neighbors(const Grid &grid)
{
  int gridm = grid.m();
  int gridn = grid.n();
  
  std::vector<int> ks, is, js, tis, tjs;
  for (int i = 0; i < gridm; ++i)
  {
    for (int j = 0; j < gridn; ++j)
    {
      int k = grid.TwoDtoOneD(i,j);
      grid.neighbors(k, ks);
      is.clear();
      js.clear();
      for (auto nk : ks)
      {
        int ni,nj;
        grid.OneDtoTwo(nk, ni,nj);
        is.push_back(ni);
        js.push_back(nj);
      }
      assert( is.size() == js.size());
      // in a periodic grid, all boxes have 20 neighbors
      assert( !grid.periodic || is.size()==20);
      
      // get neighbors another way
      tis.clear();
      tjs.clear();
      for (int jj=j-2; jj<=j+2; ++jj)
      {
        // outside grid
        auto gj = (jj + gridn) % gridn;
        if (!grid.periodic && gj != jj)
          continue;
        for (int ii=i-2; ii<=i+2; ++ii)
        {
          // outside grid
          auto gi = (ii + gridm) % gridm;
          if (!grid.periodic && gi != ii)
            continue;
          
          // center
          if (gi == i && gj == j)
            continue;
          
          // corner
          if (abs(ii-i)==2 && abs(jj-j)==2)
            continue;
          
          // good
          tis.push_back(gi);
          tjs.push_back(gj);
        }
      }
      
      // check same
      assert(tis.size()==tjs.size());
      assert(tis.size()==is.size());
      // may be in a different order, so sort to compare same
      std::set< std::pair<int,int> > tpairs, ppairs;
      for (int k = 0; k < (int) is.size(); ++k)
      {
        tpairs.emplace( tis[k], tjs[k] );
        ppairs.emplace(  is[k],  js[k] );
      }
      assert(tpairs == ppairs);
      assert(tpairs.size() == tis.size());
    } // for j
  } // for i
}

void TestGrid::test_grid(int gridm, int gridn, bool periodic, unsigned int rand_seed_p, unsigned int rand_seed_time)
{
  std::cout <<"\ntest_grid( gridm=" << gridm << ", gridn=" << gridn << ", periodic=" << periodic <<
  ", rand_seed_p=" << rand_seed_p << ", rand_seed_time=" << rand_seed_time << ")\n";

  Grid grid;
  
  grid.myrand = MyRand(rand_seed_p);
  grid.gen = std::mt19937(rand_seed_time);
  grid.periodic = periodic;
  grid.do_simple_triangulation = do_simple_triangulation;
  grid.do_early_stats = TestGrid::do_early_stats;

  grid.make_grid(gridm, gridn);
  
  grid.sweep_grid();

}


void TestGrid::test_grid_0()
{
  std::cout <<"\ntest_grid_0()\n";
  Grid grid;
  int gridm = 6;
  int gridn = 6;
  grid.periodic = false;
  grid.make_grid(gridm, gridn);
  grid.do_simple_triangulation = do_simple_triangulation;
  grid.do_early_stats = do_early_stats;
  
  TestGrid::print_all_neighbors(grid);
  TestGrid::check_all_neighbors(grid);

  grid.periodic = true;
  TestGrid::print_all_neighbors(grid);
  TestGrid::check_all_neighbors(grid);
}

void TestGrid::test_grid_2(bool periodic)
{
  std::cout <<"\ntest_grid_2()\n";
  Grid grid;
  int gridm = 6;
  int gridn = 6;
  grid.myrand = MyRand(123409234);
  grid.gen = std::mt19937(23423940);
  grid.periodic = periodic;
  grid.make_grid(gridm, gridn);
  grid.do_simple_triangulation = do_simple_triangulation;
  grid.do_early_stats = TestGrid::do_early_stats;

  // manually assign times and positions
  const auto a = sqrt(2.) / 2.;
  const auto a2 = a / 2.;
  for (int i=0; i<gridm; ++i)
  {
    for (int j=0; j<gridn; ++j)
    {
      grid.grid_square(i,j).p = Point(i*a+a2,j*a+a2);
      grid.grid_square(i,j).time = i + j + ((double) i * 0.10);
    }
  }

  std::cout << "corner-sweep manual times" << std::endl;
  grid.print_times();

  std::vector<int> early_squares;
  grid.find_early(early_squares);

  std::cout << "corner-sweep early" << std::endl;
  grid.print_early();

  std::vector<int> next_early_squares;
  next_early_squares.reserve(grid.grid_squares.size());

  std::cout << "accepting all early samples" << std::endl;
  for (int k : early_squares)
  {
    int i,j;
    grid.OneDtoTwo(k,i,j);
    assert( grid.locally_early(k) ); // expensive double-check for debugging
    grid.accept_sample(k, next_early_squares);
  }

  std::cout << "new times" << std::endl;
  grid.print_times();

  std::cout << "updated early" << std::endl;
  grid.print_early();

}

void TestGrid::test_grid_1(bool periodic)
{
  std::cout <<"\ntest_grid_1()\n";
  Grid grid;
  int gridm = 6;
  int gridn = 6;
  grid.myrand = MyRand(123409234);
  grid.gen = std::mt19937(23423940);
  grid.periodic = periodic;
  grid.do_simple_triangulation = do_simple_triangulation;
  grid.do_early_stats = TestGrid::do_early_stats;

  grid.make_grid(gridm, gridn);

  // find initially-early squares
  std::vector<int> early_squares;
  grid.find_early(early_squares);

  grid.print_times();
  grid.print_early();

}

void TestGrid::test_grid_3(bool periodic)
{
  std::cout <<"\ntest_grid_3()";
  std::cout <<" " << (periodic? "periodic" : "chopped") << "\n";
  Grid grid;
  int gridm = 6;
  int gridn = 8;
  grid.myrand = MyRand(123409234);
  grid.gen = std::mt19937(23423940);
  grid.periodic = periodic;
  grid.do_simple_triangulation = do_simple_triangulation;
  grid.do_early_stats = TestGrid::do_early_stats;

//  grid.do_print = true;
//  grid.do_counts = true;
//  grid.do_progress = true;
//  grid.do_plot = true;
  grid.do_extra_checks = true;
  // grid.do_plot_sampletime = true;
  
  grid.make_grid(gridm, gridn);
  
  // manually assign times and positions
  const auto a = sqrt(2.) / 2.;
  const auto a2 = a / 2.;
  for (int i=0; i<gridm; ++i)
  {
    for (int j=0; j<gridn; ++j)
    {
      grid.grid_square(i,j).p = Point(i*a+a2,j*a+a2);
      grid.grid_square(i,j).time = i + j + ((double) i * 0.10);
    }
  }
  
  std::cout << "corner-sweep manual times" << std::endl;
  
  grid.sweep_grid();
}

void TestGrid::test_grid_4(bool periodic)
{
  std::cout <<"\ntest_grid_4()";
  std::cout <<" " << (periodic? "periodic" : "chopped") << "\n";
  Grid grid;
  int gridm = 6;
  int gridn = 6;
  grid.myrand = MyRand(123409234);
  grid.gen = std::mt19937(23423940);
  grid.periodic = periodic;
  grid.do_simple_triangulation = do_simple_triangulation;
  grid.do_early_stats = TestGrid::do_early_stats;

//  grid.do_progress = true;
//  grid.do_counts = true;
//  // grid.do_iter_timing = true;
//  grid.do_plot = true;
//  grid.do_plot_sampletime = true;
//  grid.do_print = true;
//  grid.do_extra_checks = true;
  
  grid.make_grid(gridm, gridn);
  
  grid.sweep_grid();
}

void TestGrid::test_grid_5(bool periodic)
{
  std::cout <<"\ntest_grid_5()";
  std::cout <<" " << (periodic? "periodic" : "chopped") << "\n";
  Grid grid;
  int gridm = 60;
  int gridn = 90;
  grid.myrand = MyRand(134539234);
  grid.gen = std::mt19937(642633940);
  grid.periodic = periodic;
  grid.do_simple_triangulation = do_simple_triangulation;
  grid.do_early_stats = TestGrid::do_early_stats;

  grid.make_grid(gridm, gridn);
  
  grid.sweep_grid();
}

void TestGrid::test_grid_6(bool periodic)
{
  std::cout <<"\ntest_grid_6()";
  std::cout <<" " << (periodic? "periodic" : "chopped") << "\n";
  Grid grid;
  int gridm = 600;
  int gridn = 900;
  grid.myrand = MyRand(34539234);
  grid.gen = std::mt19937(362633940);
  grid.periodic = periodic;
  grid.do_extra_checks = false;
  grid.do_simple_triangulation = do_simple_triangulation;
  grid.do_early_stats = TestGrid::do_early_stats;

  grid.make_grid(gridm, gridn);
  
  grid.sweep_grid();
}

void TestGrid::test_grid_early(bool periodic)
{
  std::cout <<"\ntest_grid_early()";
  std::cout <<" " << (periodic? "periodic" : "chopped") << "\n";

  std::clock_t startcputime = std::clock();
  Grid grid;
  int gridm = 600;
  int gridn = 900;
  grid.myrand = MyRand(34539234);
  grid.gen = std::mt19937(362633940);
  grid.periodic = periodic;
  grid.do_simple_triangulation = true;
  grid.do_early_stats = TestGrid::do_early_stats;

  grid.do_extra_checks = false;
  grid.do_progress = true;
  grid.do_counts = true;
  grid.do_iter_timing = true;
  // grid.do_print = true;

  grid.make_grid(gridm, gridn);
  
  auto num_samples = grid.sweep_grid();

  double cpu_duration = (std::clock() - startcputime) / (double)CLOCKS_PER_SEC;

  long cells = gridm*gridn;
  printf("test_grid_early %8ld %8ld   %10.3e %7.0f %6.0f  %5.3f\n", cells, num_samples, cpu_duration, cells / cpu_duration, num_samples / cpu_duration, (double)cells/ (double)num_samples );

}

void TestGrid::test_grid_7(bool periodic)
{
  std::cout <<"\ntest_grid_7() 8.64M grid squares";
  std::cout <<" " << (periodic? "periodic" : "chopped") << "\n";
  Grid grid;
  int gridm = 600 * 4;
  int gridn = 900 * 4;
  grid.myrand = MyRand(24453417);
  grid.gen = std::mt19937(26345442);
  grid.periodic = periodic;
  grid.do_simple_triangulation = do_simple_triangulation;
  grid.do_early_stats = TestGrid::do_early_stats;

  grid.make_grid(gridm, gridn);
  
  grid.sweep_grid();
}

void TestGrid::test_grid_8(bool periodic)
{
  std::cout <<"\ntest_grid_8() multitest";
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
      Grid grid;
      unsigned seed = 24453417+ 613*gridm + 227*gridn;
      grid.myrand = MyRand(seed);
      unsigned seedtime = 26345442 + 3209*gridm + 7907*gridn;
      grid.gen = std::mt19937(seedtime);
      grid.periodic = periodic;
      grid.do_simple_triangulation = do_simple_triangulation;
      grid.do_early_stats = TestGrid::do_early_stats;

      grid.make_grid(gridm, gridn);
      grid.sweep_grid();
    }
  }
  std::cout << std::endl;
}

void TestGrid::test_grid_8_isolate1()
{
  bool periodic = false;
  std::cout <<"\ntest_grid_8_isolate";
  std::cout <<" " << (periodic? "periodic" : "chopped") << "\n";
  int gridm = 27;
  int gridn = 2;
  std::cout <<"(" << gridm << "," << gridn << ") ";
  Grid grid;
  unsigned seed = 24453417+ 613*gridm + 227*gridn;
  grid.myrand = MyRand(seed);
  unsigned seedtime = 26345442 + 3209*gridm + 7907*gridn;
  grid.gen = std::mt19937(seedtime);
  grid.periodic = periodic;
  grid.do_simple_triangulation = do_simple_triangulation;
  grid.do_early_stats = TestGrid::do_early_stats;

  grid.make_grid(gridm, gridn);
  grid.sweep_grid();
}

#include "PolygonSample.hpp"
void TestGrid::test_poly_degeneracy_1()
{
  std::cout.precision(std::numeric_limits<double>::max_digits10); // 17
  GridSquare g;

  // square
//  3.535533905932738 4.949747468305833 moveto
//  3.535533905932738 5.656854249492381 lineto
//  4.242640687119286 5.656854249492381 lineto
//  4.242640687119286 4.949747468305833 lineto
  
  const auto a = sqrt(0.5);
  const int i = round(3.535533905932738 / a); //5
  const int j = round(4.949747468305833 / a); //7
  g.polys.emplace_back( square_to_poly( (double) i*a, (double) j*a ) );

  const int gridm = 6;
  const int gridn = 8;
  const int k = i + j*gridm;
  std::cout << "Square is ersatz (" << i << "," << j << ") k=" << k << " of " <<gridm << "x" << gridn << " grid." << std::endl;

  // scoop

  // getting
//  0: circle center (3.1819805153394598,6.0104076400856501) , angles (-69.295188945364302,-44.999999999999801) (3.8890872965260099,5.3033008588991049) (3.5355339059327378,5.0749932933921666)
//  1: yminus (3.5355339059327378,5.0749932933921666) (3.5355339059327378,4.9497474683058327)
//  2: xplus (3.5355339059327378,4.9497474683058327) (3.6607797310190739,4.9497474683058327)
//  3: circle center (4.5961940777125596,4.5961940777125596) , angles (135.00000000000017,159.2951889453646) (3.6607797310190739,4.9497474683058327) (3.8890872965260099,5.3033008588991049)
//  4: circle center (4.5961940777125596,6.0104076400856501) , angles (-135.00000000000017,-135.00000000000017) (3.8890872965260099,5.3033008588991049) (3.8890872965260099,5.3033008588991049)

  // want
//  0: circle center (3.1819805153394642,6.0104076400856545) , angles (-69.295188945364586,-44.999999999999979) (3.8890872965260117,5.3033008588991075) (3.5355339059327378,5.0749932933921693)
//  1: yminus (3.5355339059327378,5.0749932933921693) (3.5355339059327378,4.9497474683058327)
//  2: xplus (3.5355339059327378,4.9497474683058327) (3.6607797310190739,4.9497474683058327)
//  3: circle center (4.5961940777125596,4.5961940777125596) , angles (135.00000000000003,159.2951889453646) (3.6607797310190739,4.9497474683058327) (3.8890872965260121,5.3033008588991066)
//  4: circle center (4.5961940777125596,6.0104076400856545) , angles (-135.00000000000006,-135) (3.8890872965260121,5.3033008588991066) (3.8890872965260117,5.3033008588991075)
  
  ScoopWorkspace workspace;
  Point c1( Point(3.1819805153394642,6.0104076400856545) );
  Point c2( Point(4.5961940777125596,4.5961940777125596) );
  Point c3( Point(4.5961940777125596,6.0104076400856545) );
  // 123
  // 132
  // 213
  // 231 // right order, wrong coordinate and angle details
  // 312 // c2 is not part of the poly
  // 321 // right order, wrong coordinate and angle details

  scoop_polys( g.polys, c3, workspace );
  scoop_polys( g.polys, c2, workspace );
  scoop_polys( g.polys, c1, workspace );
  
  struct TriangulatedPoly
  {
    double poly_area;
    Chocks chocks;
    Loop loop;
    Triangles triangles;
  };
  
  Polygon &poly = g.polys.front();
  TriangulatedPoly tp;
  auto &chocks = tp.chocks;
  auto &loop = tp.loop;
  auto &triangles = tp.triangles;
  
  // plot(poly, "square_" + std::to_string(gi) + "_" + std::to_string(gj)); // debug
  shave_chocks(poly, chocks, loop);
  std::cout << "Shaved loop\n";
  print_loop(loop);
  simplify_loop(loop);
  std::cout << "Simplified loop\n";
  print_loop(loop);

  /* int error_code = */ triangulate_loop(loop,triangles);
  tp.poly_area = polygon_area(chocks, loop, triangles);

  // no need to resample yet
  // sample_polygon(chocks, loop, triangles, gs.p, myrand, tp.poly_area);
  
  double dotsize;
  const double scale = compute_scale(gridm, gridn, dotsize);
  print_poly(poly);
  plot_g(g.p, g.time, g.polys, scale, dotsize, i,j,"test_poly_degeneracy_1");
  Chocks nochocks;
  Triangles notriangles;
  Loop noloop;
  plot_g(g.p, g.time, scale, dotsize, i,j, chocks, loop, notriangles, "test_poly_degeneracy_1-chocksonly");
  plot_g(g.p, g.time, scale, dotsize, i,j, nochocks, loop, notriangles, "test_poly_degeneracy_1-looponly");
  plot_g(g.p, g.time, scale, dotsize, i,j, nochocks, loop, triangles, "test_poly_degeneracy_1-trianglesonly");
  plot_g(g.p, g.time, scale, dotsize, i,j,chocks,loop,triangles,"test_poly_degeneracy_1");
  
}


void TestGrid::test_grid_9(bool periodic, int maxt)
{

  std::cout <<"\ntest_grid_9() linear time demo";
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

    Grid grid;
    unsigned seed = 2453417324;
    // unsigned seed = 590349834;
    grid.myrand = MyRand(seed);
    unsigned seedtime = 2634524345;
    // unsigned seedtime = 23967113;
    grid.gen = std::mt19937(seedtime);
    grid.periodic = periodic;
    grid.do_counts = false;
    grid.do_extra_checks = false;
    grid.do_simple_triangulation = do_simple_triangulation;
    grid.do_early_stats = TestGrid::do_early_stats;

    // debug
    // grid.do_counts = true;
    // grid.do_extra_checks = true;
    // grid.do_plot = true;
    // grid.do_progress = true;
    // grid.do_print = true;
    
    grid.make_grid(gridm, gridn);
    auto num_samples = grid.sweep_grid();
    
    auto endcputime = std::clock();
    // auto chrono_end = std::chrono::high_resolution_clock::now();
    double cpu_duration = (endcputime - startcputime) / (double)CLOCKS_PER_SEC;
    // double chrono_duration = std::chrono::duration_cast<std::chrono::duration<double> >(chrono_end - chrono_start).count();
    
    printf("%2d %8ld %8ld   %10.3e %7.0f %6.0f  %5.3f -clock\n", t, ag[t], num_samples, cpu_duration, ag[t] / cpu_duration, num_samples / cpu_duration, (double)ag[t]/ (double)num_samples );
    // printf("%2d %8ld %8ld   %10.3e %7.0f %6.0f  %5.3f -chrono\n", t, ag[t], num_samples, chrono_duration, ag[t] / chrono_duration, num_samples / chrono_duration, (double)ag[t]/ (double)num_samples );
  }
  std::cout << std::endl;
}

void grid_figure()
{
  bool periodic = true;
  std::cout <<"\ngrid_figure()";
  std::cout <<" " << (periodic? "periodic" : "chopped") << "\n";
  Grid grid;
  int gridm = 9;
  int gridn = 6;
  grid.myrand = MyRand(134539234);
  grid.gen = std::mt19937(62333940);
  grid.periodic = periodic;
  grid.do_simple_triangulation = false;
  grid.point_prepasses=0;
  grid.square_prepasses=0;

  grid.do_plot = true;
  grid.do_print = true;
  grid.do_extra_checks = true;
  grid.do_progress = true;
  grid.do_counts = true;
  grid.plot_dotsize_factor = 3.;
  grid.plot_linewidth_factor = 20.;

  grid.make_grid(gridm, gridn);
  grid.sweep_grid();
  
}

void DeterministicMPS_runtime_table()
{
  TestGrid::do_simple_triangulation = true;
  TestGrid::do_early_stats = false;
  TestGrid::test_grid_9(true, 16); //16
}

void PSA_point_set_analysis_figure()
{
  std::cout <<"\nPSA_point_set_analysis_figure()\n";
  std::cout << "This generates distributions in the input file format for PSA: http://code.google.com/p/psa/ \n";
  std::cout << "The domain is scaled to the box [0,1]x[0,1] on output.\n";
  std::cout << "To make the paper figures, run this method to create files `pointsK.txt' where K=[1,100].\n" <<
  "Then make psa, and do two psa runs:\n"
  // "% psa points.txt" << std::endl;
  "% psa --avg points*.txt\n" <<
  "% psa --avg --ani points*.txt" << std::endl;
  std::cout << "The first produces avg.pdf, and the second produces avg_ani.tex, which when run through latex produces avg_ani.pdf\n";
  std::cout << "Note PSA's averaging is limited to the same number of points per file, so each file is truncated to the minimum number of points over all files.\n"
  "psa will report 'Analyzing only the first 3449 points from each file' but any defects in the figures are unnoticable\n"
  "because we shuffled our points to random order instead of grid-scan order before outputting them, so which ones were truncated is spatially random.\n"
  "Also, because of scaling, in avg.pdf the reported ``mindist'' is not the true value, which is 1 up to machine precision.\n"
  << std::endl;

  //  @Misc{      psa,
  //    author  = {Thomas Schl\"omer},
  //    title    = {{PSA} Point Set Analysis},
  //    howpublished  = {Version 1.1, \url{http://code.google.com/p/psa/}},
  //    year    = {2013}
  //  }
  
  for (int i = 1; i<= 100; ++i)
  {
    Grid grid;
    int gridm = 100;
    int gridn = gridm;
    // makes about 3.5k points, PSA takes 1 minute with its O(n^2) runtime.
    
    unsigned seed = 24453417 + 613239*i;
    grid.myrand = MyRand(seed);
    unsigned seedtime = 26345442 + 89345209*i;
    grid.gen = std::mt19937(seedtime);
    
    grid.periodic = true;
    grid.do_simple_triangulation = true;
    grid.do_early_stats = false; // they'd overwrite each other

    grid.do_plot = false;
    grid.do_print = false;
    grid.do_extra_checks = true;
    grid.do_progress = false;
    grid.do_counts = true;
    grid.do_PSA_points = i;
    
    grid.make_grid(gridm, gridn);
    grid.sweep_grid();
  }
}


void DeterministicMPS( int gridm, int gridn, bool periodic )
{
  Grid grid;
  grid.periodic = periodic;

  unsigned seed = 1127308102;
  unsigned seedtime = 690798395;
  grid.myrand = MyRand(seed);
  grid.gen = std::mt19937(seedtime);
  
  grid.do_simple_triangulation = true;
  grid.do_early_stats = false;

  grid.do_plot = true;
  grid.do_print = false;
  grid.do_extra_checks = true;
  grid.do_progress = true;
  grid.do_counts = true;
  grid.do_PSA_points = false;

  grid.point_prepasses=0; // 7 is best for running, but 0 is better to illustrate the main algorithm in the paper.
  grid.square_prepasses=0;

  grid.make_grid(gridm, gridn);
  grid.sweep_grid();
}
