//
//  TestPolygonSample.cpp
//  AvoidReject
//
//  Created by Mitchell, Scott A on 10/21/21.
//  Copyright © 2021 Mitchell, Scott A. All rights reserved.
//

#include "TestPolygonSample.hpp"

#include <iostream>
#include <fstream>

#include "Geometry.hpp"
#include "Triangulate.hpp"
#include "PolygonSample.hpp"
#include "TestTriangulate.hpp"


void plot_samples(const std::vector<Point> &samples, std::ofstream &fout)
{
  fout << "% " << samples.size() << " sample points in polygon\n";
  for (size_t i=0; i<samples.size(); ++i)
  {
    const Point &p = samples[i];
    fout << "\n% sample point " << i << "\n";
    fout << "newpath " << p.x << " " << p.y << " 0.001 0 360 arc closepath fill\n";
  }
}


void test_polygon_sample(const std::vector<Polygon> &polys, MyRand &myrand, size_t num_samples, int id, int sub_id = -1, bool do_plot = true)
{
  
  Chocks chocks;
  Loop loop;
  Triangles triangles;
  for (auto &poly : polys)
  {
    shave_chocks(poly, chocks, loop);
    triangulate_loop(loop,triangles);
  
    const double poly_area = polygon_area(chocks, loop, triangles);

    // constant density if num_samples not specified
    if (num_samples == (size_t)-1)
    {
      num_samples = std::max(1.,20000. * poly_area);
    }
    std::vector<Point> samples;
    for (size_t i = 0; i < num_samples; ++i)
    {
      
      // return random point p from polygon
      //   prior to calling, shave chocks and triangulate loop
      // poly_area can be provided as a speedup, otherwise we'll calculate it from scratch
      Point p;
      sample_polygon(chocks, loop, triangles, p, myrand, poly_area);
      
      samples.push_back(p);
    } // for samples

    // plot chocks and triangles
    std::string fname = "test_polygon_sample_" + std::to_string(id);
    if (sub_id >= 0)
    {
      fname += "_" + std::to_string(sub_id);
    }
    if (do_plot)
    {
      const bool do_labels = false;
      std::ofstream fout;
      fout.open (fname + "-chocks_triangles.ps");
      plot_preamble(fout);
      plot_triangles(loop, triangles, fout, do_labels);
      plot_chocks(chocks,fout,false);
      fout << "\nshowpage\n";
      fout.close();
    }
    
    // plot samples
    if (do_plot)
    {
      std::ofstream fout;
      fout.open (fname + "-samples.ps");
      plot_preamble(fout);
      // void plot_bounding_box(Point bbox_lo, Point bbox_hi, bool expand_box, std::ofstream &fout);
      plot_samples(samples, fout);
      fout << "\nshowpage\n";
      fout.close();
    }
    
    // plot together
    if (do_plot)
    {
      const bool do_labels = false;
      std::ofstream fout;
      fout.open (fname + "-chocks_triangles_samples.ps");
      plot_preamble(fout);
      plot_triangles(loop, triangles, fout, do_labels);
      plot_chocks(chocks,fout,false);
      plot_samples(samples, fout);
      fout << "\nshowpage\n";
      fout.close();
    }

  } // for polys
  
}

void test_polygon_sample_bites(const std::vector<Point> &bites, MyRand &myrand, const size_t num_samples, int id, int sub_id = -1)
{
  
  Polygon square_poly = square_to_poly(0., 0.);
  
  std::vector<Polygon> polys = {square_poly};
  std::vector<Polygon> workspace;
  for (auto &d : bites)
  {
    scoop_polys(polys, d, workspace);
  }
  
  test_polygon_sample(polys, myrand, num_samples, id, sub_id);
}

void test_polygon_sample_0()
{
  // make a polygon
  std::vector<Point> bites =
  {
    Point(-0.4, -0.50),
    Point(0.3,   1.35),
    Point(1.05,  -0.5)
  };
  MyRand myrand(483738);
  test_polygon_sample_bites(bites, myrand, 1200, 0);
}

void test_polygon_sample_1()
{
  // make a polygon
  std::vector<Point> bites =
  {
    Point(-0.4, -0.50),
    Point(0.3,   1.35),
    Point(1.05,  -0.5)
  };
  MyRand myrand(234233738);
  test_polygon_sample_bites(bites, myrand, 1200, 1);
}
void test_polygon_sample_2()
{
  // make a polygon
  std::vector<Point> bites =
  {
    Point(-0.4, -0.50),
    Point(0.3,   1.35),
    Point(1.05,  -0.5)
  };
  MyRand myrand(423834533);
  test_polygon_sample_bites(bites, myrand, 1200, 2);
}

void test_polygon_sample_3()
{
  // sample from randomly bitten polygons as in "TestTriangulate.cpp"
  MyRand myrand(3829572716);
  const size_t num_samples = (size_t) -1; // 10000;
  
  std::vector<Polygon> polys;
  std::vector<Point> bites;
  
  for (int i=0; i<100; ++i)
  {
    polys.clear();
    bites.clear();
    int max_bites = 1+std::floor( myrand.u()*7 );
    random_bites(polys, bites, max_bites, myrand);
    test_polygon_sample(polys, myrand, num_samples, 3, i);
  }
}


void test_polygon_sample_4()
{
  // sample from randomly bitten polygons as in "TestTriangulate.cpp"
  MyRand myrand(3459576234);
  const size_t num_samples = (size_t) -1; // 10000;
  
  std::vector<Polygon> polys;
  std::vector<Point> bites;
  
  for (int i=0; i<1000000; ++i)
  {
    polys.clear();
    bites.clear();
    int max_bites = 1+std::floor( myrand.u()*7 );
    random_bites(polys, bites, max_bites, myrand);
    test_polygon_sample(polys, myrand, num_samples, 4, i, false);
  }
}
