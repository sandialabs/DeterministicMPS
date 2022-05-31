//
//  TestTriangulate.cpp
//  AvoidReject
//
//  Created by Mitchell, Scott A on 10/20/21.
//  Copyright © 2021 Mitchell, Scott A. All rights reserved.
//

#include "TestTriangulate.hpp"
#include "TestGeometry.hpp"
#include "MyRand.hpp"
#include "Grid.hpp"

#include <ostream>
#include <iostream>
#include <fstream>

void test_shave(const Polygon &poly, Chocks &chocks, Loop &loop, Point bbox_lo, Point bbox_hi, std::string fname, bool do_plots)
{
  if (do_plots)
  {
    std::string fname_poly = fname + "_poly";
    std::cout << "test_shave writing polygon to " << fname_poly << std::endl;
    plot_poly(poly, fname_poly, bbox_lo, bbox_hi);
  }
  
  shave_chocks(poly, chocks, loop);
  
  // plot chocks
  if (do_plots)
  {
    std::string fname_chock = fname + "_chocks.ps";
    std::cout << "test_shave writing chocks to " << fname_chock << std::endl;
    std::ofstream fout;
    fout.open (fname_chock);
    plot_preamble(fout);
    plot_bounding_box(bbox_lo, bbox_hi, true, fout);
    fout << "\n% " << chocks.size() << " chocks\n";
    size_t ci=0;
    for (auto &chock : chocks)
    {
      fout << "\n% chock " << ci++ << "\n";
      plot_chock(chock.c, chock.t, chock.q, fout);
    }
  }
  
  // plot loop
  if (do_plots)
  {
    std::string fname_loop = fname + "_loop.ps";
    std::cout << "test_shave writing loop to " << fname_loop << std::endl;
    std::ofstream fout;
    fout.open(fname_loop);
    plot_preamble(fout);
    plot_bounding_box(bbox_lo, bbox_hi, true, fout);
    plot_loop(loop, fout);
  }
  
}

void scoop_square(std::vector<Polygon> &polys, Point &bbox_lo, Point &bbox_hi, const std::vector<Point> &disks)
{
  double x0, y0, xa, ya, a;
  Polygon square;
  default_square_poly(x0, y0, xa, ya, a, square );
  bbox_lo = Point(x0,y0);
  bbox_hi = Point(xa,ya);

  // scoop a poly
  polys = {square};
  std::vector<Polygon> workspace;
  for (auto &d : disks)
  {
    scoop_polys(polys, d, workspace);
  }
}

void plot_for_paper(const Chocks &chocks, Loop &loop, const Triangles &triangles, std::string fname, double linewidth = 0.005)
{
  double scale, dotsize;
  scale = compute_scale(3, 3, dotsize);
  
  fname += ".ps";
  std::cout << "plotting " << fname  << std::endl;
  std::ofstream fout;
  fout.open(fname);
    
  // scale
  grid_preamble(scale, fout);
  // bigger line for figure clarity
  fout << linewidth << " setlinewidth\n";


  Point bbox_lo(0,0);
  Point bbox_hi( sqrt(0.5), sqrt(0.5) );
  plot_bounding_box(bbox_lo, bbox_hi, false, fout);
  plot_loop(loop, fout, false, true );
  plot_chocks(chocks, fout);
  plot_full_circles(chocks, fout);
  plot_triangles(loop, triangles, fout, false);
  
  fout << "\nshowpage\n";
  fout.close();
}


void test_shave(const std::vector<Point> &bites, bool do_triangulate, size_t id, size_t id2 = -1, bool do_plots = true)
{
  std::vector<Polygon> polys;
  Point bbox_lo, bbox_hi;
  scoop_square(polys, bbox_lo, bbox_hi, bites);
  
  std::string fname = "test_shave_" + std::to_string(id);
  if (id2 != -1)
  {
    fname += "_sub" + std::to_string(id2);
  }

  if (do_plots)
  {
    plot_square_biters( bbox_lo, bbox_hi, bites, bites.size()-1, fname );
  }
  
  Chocks chocks;
  Loop loop;
  for (size_t pi = 0; pi < polys.size(); ++pi)
  {
    std::string namestring = fname + "_" + std::to_string(pi);    
    test_shave(polys[pi], chocks, loop, bbox_lo, bbox_hi, namestring, do_plots);
    
    if (do_triangulate)
    {
      Triangles triangles;
      // triangulate_loop_simple(loop,triangles); // 0.6s
      triangulate_loop(loop,triangles); // 4.4s

      if (do_plots)
      {
        plot_triangulation(loop, triangles, true, bbox_lo, bbox_hi, namestring + "_triangles");
        
        
        // plot bites and chocks and triangulation together
      }
    }
    chocks.clear();
    loop.clear();
  }
}

void test_st_0(bool do_triangulate)
{
  std::vector<Point> bites =
  {
    Point(-0.4, -0.60),
    Point(0.2,   1.2),
    Point(0.7,  -0.7)
  };
  test_shave(bites,do_triangulate,0);
}

void test_st_1(bool do_triangulate)
{
  std::vector<Point> bites =
  {
    Point(-0.4, -0.50),
    Point(0.3,   1.35),
    Point(1.05,  -0.5)
  };
  test_shave(bites,do_triangulate,1);
}

void test_st_2(bool do_triangulate)
{
  MyRand myrand(3829572716);

  std::vector<Polygon> polys;
  std::vector<Point> bites;
  
  for (size_t i=0; i<100; ++i)
  {
    if (i==95)
    {
      std::cout << "debug me" << std::endl;
    }
    polys.clear();
    bites.clear();
    int max_bites = 1+std::floor( myrand.u()*7 );
    random_bites(polys, bites, max_bites, myrand);
    test_shave(bites,do_triangulate,2,i);
  }
}

void test_st_3(bool do_triangulate)
{
  // time: 1,000,000 points takes about
  //   34 seconds with fancy triangulation, compiled optimized
  //   simple ear cutting
  MyRand myrand(234938716);
  
  std::vector<Polygon> polys;
  std::vector<Point> bites;
  
  size_t maxi = 1000000;
  for (size_t i=0; i<maxi; ++i)
  {
    polys.clear();
    bites.clear();
    int max_bites = 1+std::floor( myrand.u()*7 );
    random_bites(polys, bites, max_bites, myrand);
    test_shave(bites,do_triangulate,3,i,false);
  }
}

void test_shave_0()
{
  test_st_0(false);
}
void test_shave_1()
{
  test_st_1(false);
}
void test_shave_2()
{
  test_st_2(false);
}
void test_shave_3()
{
  test_st_3(false);
}


void test_triangulate_0()
{
  test_st_0(true);
}
void test_triangulate_1()
{
  test_st_1(true);
}
void test_triangulate_2()
{
  test_st_2(true);
}
void test_triangulate_3()
{
  test_st_3(true);
}


void test_angle_measure()
{
  Point a(1,0), b(0,0), c(1,0);
  std::cout << "  " << "theta" << " \t" << "angle_measure" << "\n";
  for (double theta = 0; theta<=2*M_PI; theta+= M_PI * 0.01)
  {
    a.x = cos(theta);
    a.y = sin(theta);
    auto ms = measure_angle(a,b,c);
    std::cout << "  " << theta << " \t" << ms << "\n";
  }
  std::cout << "expected to be monotonically increasing " << std::endl;
}

void test_circum()
{
  Point a(0,1), b(0,0), c(1,0);
  Point o;
  double R;
  circumcenter(a, b, c, o, R);
  Point m(0.6,0.6);
  auto is_empty = is_empty_circle_b(a, b, c, m, o, R);
  std::cout << "\ncircumcenter: (" << o.x << ", " << o.y << ") radius: " << R << std::endl;
  std::cout << "test point m: (" << m.x << ", " << m.y << ") empty_circle_b = " << is_empty << std::endl;
  
  m=Point(-0.1,-0.1);
  is_empty = is_empty_circle_b(a, b, c, m, o, R);
  std::cout << "\ncircumcenter: (" << o.x << ", " << o.y << ") radius: " << R << std::endl;
  std::cout << "test point m: (" << m.x << ", " << m.y << ") empty_circle_b = " << is_empty << std::endl;
  
  m=Point(-0.01,0.3);
  is_empty = is_empty_circle_b(a, b, c, m, o, R);
  std::cout << "\ncircumcenter: (" << o.x << ", " << o.y << ") radius: " << R << std::endl;
  std::cout << "test point m: (" << m.x << ", " << m.y << ") empty_circle_b = " << is_empty << std::endl;
  
  m=Point(0.01,0.3);
  is_empty = is_empty_circle_b(a, b, c, m, o, R);
  std::cout << "\ncircumcenter: (" << o.x << ", " << o.y << ") radius: " << R << std::endl;
  std::cout << "test point m: (" << m.x << ", " << m.y << ") empty_circle_b = " << is_empty << std::endl;
}

void test_grazing_0()
{
  // same y coordinate
  const double a = sqrt(2.)/2.; // side length
  std::vector<Point> bites =
  {
    Point(-0.5, a/2. /* + std::numeric_limits<double>::epsilon() */),
    Point(1.6,  a/2.)
  };
  // approach kissing
  for (auto i = 0; i<100; ++i)
  {
    bites[1].x = 1.5 + 1./std::pow(2, i);
    test_shave(bites,true, 5, i, true);
  }
  // eps separation to eps crossing
  for (auto i = 0; i<8; ++i)
  {
    bites[1].x = 1.5 + ((3-i)*std::numeric_limits<double>::epsilon());
    test_shave(bites,true, 6, i, true);
  }
}

void test_grazing_1()
{
  
  const double a = sqrt(2.)/2.; // side length
  std::vector<Point> bites =
  {
    Point(-0.3, -0.3),
    Point(1.6,  a/2.)
  };
  double theta = (M_PI / 2.) * 0.32;
  for ( ; theta < (M_PI / 2.); theta += (M_PI / 2.) * 0.01)
  {
    // approach kissing
    for (auto i = 0; i<100; ++i)
    {
      double center_dist = 2. + 1./std::pow(2, i);
      bites[1].x = bites[0].x + center_dist * std::cos(theta);
      bites[1].y = bites[0].y + center_dist * std::sin(theta);
      test_shave(bites,true, 7, i, false);
    }
    // eps separation to eps crossing
    for (auto i = 0; i<8; ++i)
    {
      double center_dist = 2. + ((3-i)*std::numeric_limits<double>::epsilon());
      bites[1].x = bites[0].x + center_dist * std::cos(theta);
      bites[1].y = bites[0].y + center_dist * std::sin(theta);
      test_shave(bites,true, 8, i, false);
    }
  }
}

void test_grazing_2()
{
  // same y coordinate
  const double a = sqrt(2.)/2.; // side length
  std::vector<Point> bites =
  {
    Point(1., a/2. /* + std::numeric_limits<double>::epsilon() */)
  };
  // approach kissing
  for (auto i = 0; i<100; ++i)
  {
    bites[0].x = 1. + 1./std::pow(2, i);
    test_shave(bites, true, 9, i, true);
  }
  // eps separation to eps crossing
  for (auto i = 0; i<8; ++i)
  {
    const auto delta = (3-i)*std::numeric_limits<double>::epsilon();
    bites[0].x = 1. + delta;
    test_shave(bites,true, 10, i, true);
  }
}

void figure_trim(std::vector<Point> &bites, std::string fname, bool do_triangles, double linewidth = 0.005)
{
  std::vector<Polygon> polys;
  Point bbox_lo, bbox_hi;
  scoop_square(polys, bbox_lo, bbox_hi, bites);
  
  Chocks chocks;
  Loop loop;
  Triangles triangles, no_triangles;

  for (size_t pi = 0; pi < polys.size(); ++pi)
  {
    Polygon poly = polys[pi];
    
    shave_chocks(poly, chocks, loop);

    // triangulate_loop_simple(loop,triangles);
    triangulate_loop(loop,triangles);
        
    plot_for_paper(chocks, loop,
                   (do_triangles ? triangles : no_triangles),
                   fname + (do_triangles ? "_triangles" : "_" ) + std::to_string(pi),
                   linewidth);
    
    chocks.clear();
    loop.clear();
  }
}

void figure_trimA(bool do_triangles = false)
{
  std::vector<Point> bites =
  {
    Point(-0.4, -0.2),
    Point(1.5,  -0.3)
  };
  figure_trim(bites, "fig_trim_A", do_triangles);
}

void figure_trimB(bool do_triangles = false)
{
  std::vector<Point> bites =
  {
    Point(0.3, -0.35)
  };
  figure_trim(bites, "fig_trim_B", do_triangles);
}

void figure_trimC(bool do_triangles = false)
{
  std::vector<Point> bites =
  {
    Point(-0.4, -0.4),
    Point(1.1, 1.1)
  };
  figure_trim(bites, "fig_trim_C", do_triangles);
}

void figure_trimD(bool do_triangles = false)
{
  std::vector<Point> bites =
  {
    Point(-0.658, 0.35),
    Point(0.75, -0.85),
    Point(1.358, 0.35)
  };
  figure_trim(bites, "fig_trim_D", do_triangles, 0.001);
}

void figure_trimE(bool do_triangles = false)
{
  std::vector<Point> bites =
  {
    Point(-0.644, 0.35),
    Point(0.75, -0.85),
    Point(1.350, 0.35)
  };
  figure_trim(bites, "fig_trim_E", do_triangles, 0.001);
}


void figure_trim()
{
  figure_trimA();
  figure_trimB();
  figure_trimC();
}

void figure_triangles()
{
  figure_trimA(true);
  figure_trimB(true);
  figure_trimC(true);
  figure_trimD(true);
  figure_trimE(true);
}
