//
//  TestGeometry.cpp
//  AvoidReject
//
//  Created by Mitchell, Scott A on 8/20/21.
//  Copyright © 2021 Mitchell, Scott A. All rights reserved.
//

#include "TestGeometry.hpp"
#include <assert.h>
#include <string>

#include "Geometry.hpp"
#include "MyRand.hpp"

void test_intersect_circle()
{
  
  bool disjoint;
  Point n0, n1;

  Point d0( Point(0.,0.) );
  Point d1( Point(1.,0.) );

  intersect_circle(d0, d1, disjoint, n0, n1);
  assert(!disjoint);
  assert( fabs(n0.x - n1.x) < 10 * __DBL_EPSILON__);
  assert( fabs(n0.x - 0.5 ) < 1000 * __DBL_EPSILON__);
  assert( fabs( n0.y + n1.y ) < 10 * __DBL_EPSILON__);
  assert( fabs(n0.y) - sqrt(3)/2. < 1000 * __DBL_EPSILON__);
  assert( fabs(n1.y) - sqrt(3)/2. < 1000 * __DBL_EPSILON__);

  d1.x = -1.;
  d1.y = 0.;
  intersect_circle(d0, d1, disjoint, n0, n1);
  assert(!disjoint);
  assert( fabs(n0.x - n1.x) < 10 * __DBL_EPSILON__);
  assert( fabs(n0.x + 0.5 ) < 1000 * __DBL_EPSILON__);
  assert( fabs( n0.y + n1.y ) < 10 * __DBL_EPSILON__);
  assert( fabs(n0.y) - sqrt(3)/2. < 1000 * __DBL_EPSILON__);
  assert( fabs(n1.y) - sqrt(3)/2. < 1000 * __DBL_EPSILON__);

  d1.x = 2.0;
  d1.y = 0.6;
  intersect_circle(d0, d1, disjoint, n0, n1);
  assert(disjoint);

  d1.x = 1.9;
  d1.y = 0.6;
  intersect_circle(d0, d1, disjoint, n0, n1);
  assert(!disjoint);

  d1.x = 0.;
  d1.y = 1.;
  intersect_circle(d0, d1, disjoint, n0, n1);
  assert(!disjoint);
  assert( fabs(n0.y - n1.y) < 10 * __DBL_EPSILON__);
  assert( fabs(n0.y - 0.5 ) < 1000 * __DBL_EPSILON__);
  assert( fabs( n0.x + n1.x ) < 10 * __DBL_EPSILON__);
  assert( fabs(n0.x) - sqrt(3)/2. < 1000 * __DBL_EPSILON__);
  assert( fabs(n1.x) - sqrt(3)/2. < 1000 * __DBL_EPSILON__);
  
  d1.x = 0.;
  d1.y = -1.;
  intersect_circle(d0, d1, disjoint, n0, n1);
  assert(!disjoint);
  assert( fabs(n0.y - n1.y) < 10 * __DBL_EPSILON__);
  assert( fabs(n0.y + 0.5 ) < 1000 * __DBL_EPSILON__);
  assert( fabs( n0.x + n1.x ) < 10 * __DBL_EPSILON__);
  assert( fabs(n0.x) - sqrt(3)/2. < 1000 * __DBL_EPSILON__);
  assert( fabs(n1.x) - sqrt(3)/2. < 1000 * __DBL_EPSILON__);
  
  // tests from drawings in powerpoint, 4 different quadrants
  Point N0, N1; // expected answers

  d0 = Point(6.91,5.04);
  d1 = Point(7.61,6.18);
  N0 = Point(7.89,5.21);
  N1 = Point(6.62,5.99);
  intersect_circle(d0, d1, disjoint, n0, n1);
  if (n0.y>n1.y)
    std::swap(n0,n1);

  assert(!disjoint);
  assert( fabs(n0.x - N0.x) < 0.02 );
  assert( fabs(n0.y - N0.y) < 0.02 );
  assert( fabs(n1.x - N1.x) < 0.02 );
  assert( fabs(n1.y - N1.y) < 0.02 );

  //
  d0 = Point(1.6,2.35);
  d1 = Point(2.98,1.48);
  N0 = Point(1.98,1.41);
  N1 = Point(2.6,2.4);
  intersect_circle(d0, d1, disjoint, n0, n1);
  if (n0.y>n1.y)
    std::swap(n0,n1);
  
  assert(!disjoint);
  assert( fabs(n0.x - N0.x) < 0.02 );
  assert( fabs(n0.y - N0.y) < 0.02 );
  assert( fabs(n1.x - N1.x) < 0.02 );
  assert( fabs(n1.y - N1.y) < 0.02 );
  
  //
  d0 = Point(7.06,2.53);
  d1 = Point(6.13,1.72);
  N0 = Point(7.1,1.53);
  N1 = Point(6.08,2.72);
  intersect_circle(d0, d1, disjoint, n0, n1);
  if (n0.y>n1.y)
    std::swap(n0,n1);
  
  assert(!disjoint);
  assert( fabs(n0.x - N0.x) < 0.02 );
  assert( fabs(n0.y - N0.y) < 0.02 );
  assert( fabs(n1.x - N1.x) < 0.02 );
  assert( fabs(n1.y - N1.y) < 0.02 );
  
  //
  d0 = Point(3.64,4.91);
  d1 = Point(2.01,5.89);
  N0 = Point(2.66,5.13);
  N1 = Point(2.98,5.66);
  intersect_circle(d0, d1, disjoint, n0, n1);
  if (n0.y>n1.y)
    std::swap(n0,n1);
  
  assert(!disjoint);
  assert( fabs(n0.x - N0.x) < 0.02 );
  assert( fabs(n0.y - N0.y) < 0.02 );
  assert( fabs(n1.x - N1.x) < 0.02 );
  assert( fabs(n1.y - N1.y) < 0.02 );
  
}

void test_intersect_x()
{
  bool disjoint;
  Point n0, n1;
 
  Circle d0( Point(0.,0.) );
  intersect_x(d0, 0, disjoint, n0, n1);
  assert( fabs(n0.x) < 10 * __DBL_EPSILON__);
  assert( fabs(n1.x) < 10 * __DBL_EPSILON__);
  assert( fabs(n0.y + 1) < 10 * __DBL_EPSILON__);
  assert( fabs(n1.y - 1) < 10 * __DBL_EPSILON__);
  assert(!disjoint);
  
  intersect_x(d0, -1.3, disjoint, n0, n1);
  assert(disjoint);
  intersect_x(d0, 1.01, disjoint, n0, n1);
  assert(disjoint);

  d0.c = Point(12.,-4.);
  intersect_x(d0, 11.3, disjoint, n0, n1);
  assert( fabs(n0.x - 11.3) < 10 * __DBL_EPSILON__);
  assert( fabs(n1.x - 11.3) < 10 * __DBL_EPSILON__);
  const double expect_y0 = -4 - sqrt( 1 - 0.7*0.7 );
  const double expect_y1 = -4 + sqrt( 1 - 0.7*0.7 );
  assert( fabs(n0.y - expect_y0 ) < 100 * __DBL_EPSILON__);
  assert( fabs(n1.y - expect_y1 ) < 100 * __DBL_EPSILON__);

}
void test_intersect_y()
{
  bool disjoint;
  Point n0, n1;
  
  Circle d0( Point(0.,0.) );
  intersect_y(d0, 0, disjoint, n0, n1);
  assert( fabs(n0.y) < 10 * __DBL_EPSILON__);
  assert( fabs(n1.y) < 10 * __DBL_EPSILON__);
  assert( fabs(n0.x + 1) < 10 * __DBL_EPSILON__);
  assert( fabs(n1.x - 1) < 10 * __DBL_EPSILON__);
  assert(!disjoint);

  d0.c = Point(3.,1.);
  intersect_y(d0, 1.2, disjoint, n0, n1);
  assert( fabs(n0.y - 1.2) < 10 * __DBL_EPSILON__);
  assert( fabs(n1.y - 1.2) < 10 * __DBL_EPSILON__);
  const double expect_x0 = 3 - sqrt( 1 - 0.2*0.2 );
  const double expect_x1 = 3 + sqrt( 1 - 0.2*0.2 );
  assert( fabs(n0.x - expect_x0) < 100 * __DBL_EPSILON__);
  assert( fabs(n1.x - expect_x1) < 100 * __DBL_EPSILON__);
  assert(!disjoint);

  intersect_y(d0, -0.1, disjoint, n0, n1);
  assert(disjoint);

}

void test_polygon0()
{
  
  const auto a = 1 / sqrt(2);
  Polygon poly =
  {
    Segment(Point(0,0), Point(a,0), XPLUS, Point(0,0)),
    Segment(Point(a,0), Point(a,a), YPLUS, Point(a,0)),
    Segment(Point(a,a), Point(0,a), XMINUS, Point(a,a)),
    Segment(Point(0,a), Point(0,0), YMINUS, Point(0,a))
  };
  // plot_poly( poly, "unitsquare" );
  
  if (/* DISABLES CODE */ (1))
  {
    std::vector<Polygon> result;
    scoop_poly(poly, Point(-0.3,-0.50), result);
    size_t i = 0;
    for (auto &r : result)
    {
      plot_poly( r, "cornerbite-one" + std::to_string(i) );
      ++i;
    }
  }

  if (/* DISABLES CODE */ (1))
  {
    std::vector<Polygon> result;
    scoop_poly(poly, Point(0.5/sqrt(2.),-0.270), result);
    size_t i = 0;
    for (auto &r : result)
    {
      plot_poly( r, "cornerbite-zero" + std::to_string(i) );
      ++i;
    }
  }

  if (/* DISABLES CODE */ (1))
  {
    std::vector<Polygon> result;
    scoop_poly(poly,  Point(0.5/sqrt(2.),-0.98), result);
    size_t i = 0;
    for (auto &r : result)
    {
      plot_poly( r, "cornerbite-zerolow" + std::to_string(i) );
      ++i;
    }
  }

  if (/* DISABLES CODE */ (1))
  {
    std::vector<Polygon> result;
    scoop_poly(poly, Point(0.1,0.9), result);
    size_t i = 0;
    for (auto &r : result)
    {
      plot_poly( r, "cornerbite-three" + std::to_string(i) );
      ++i;
    }
  }

  if (/* DISABLES CODE */ (1))
  {
    std::vector<Polygon> result;
    scoop_poly(poly, Point(0.3,0.4), result);
    size_t i = 0;
    for (auto &r : result)
    {
      plot_poly( r, "cornerbite-four" + std::to_string(i) );
      ++i;
    }
  }

  if (/* DISABLES CODE */ (1))
  {
    std::vector<Polygon> result;
    scoop_poly(poly, Point(2.,0.0), result);
    size_t i = 0;
    for (auto &r : result)
    {
      plot_poly( r, "cornerbite-miss" + std::to_string(i) );
      ++i;
    }
  }

}

void test_polygon1()
{
  const double x0 = 1.2;
  const double y0 = -0.4;
  const auto a = 1 / sqrt(2);
  const double xa = x0+a;
  const double ya = y0+a;
  const Polygon poly =
  {
    Segment(Point(x0,y0), Point(xa,y0), XPLUS,  Point(x0,y0)),
    Segment(Point(xa,y0), Point(xa,ya), YPLUS,  Point(xa,y0)),
    Segment(Point(xa,ya), Point(x0,ya), XMINUS, Point(xa,ya)),
    Segment(Point(x0,ya), Point(x0,y0), YMINUS, Point(x0,ya))
  };

  const double xcmin = x0 - 1.1;
  const double xcmax = xa + 1.1;
  const double ycmin = y0 - 1.1;
  const double ycmax = ya + 1.1;
  
  const double dx = 0.111;
  const double dy = 0.092;
  size_t j = 0;
  for (double x=xcmin; x<=xcmax; x+=dx)
  {
    for (double y=ycmin; y<=ycmax; y+=dy)
    {
      Polygon polycopy = poly;
      std::vector<Polygon> result;
      scoop_poly(polycopy, Point(x,y), result);
      
      size_t i = 0;
      for (auto &r : result)
      {
        plot_poly( r, "test_polygon1_" + std::to_string(j) + "_" + std::to_string(i) );
        ++i;
      }
      ++j;
    }
  }
}

void test_polygon2()
{
  const double x0 = -0.2;
  const double y0 = -0.3;
  const auto a = 1 / sqrt(2);
  const double xa = x0+a;
  const double ya = y0+a;
  const Polygon poly =
  {
    Segment(Point(x0,y0), Point(xa,y0), XPLUS,  Point(x0,y0)),
    Segment(Point(xa,y0), Point(xa,ya), YPLUS,  Point(xa,y0)),
    Segment(Point(xa,ya), Point(x0,ya), XMINUS, Point(xa,ya)),
    Segment(Point(x0,ya), Point(x0,y0), YMINUS, Point(x0,ya))
  };

  MyRand myrand(234903240);
  
  const size_t num_tests = 1000;
  for (size_t j = 0; j < num_tests; ++j)
  {
    auto randx = 1.5-3.*myrand.u();
    auto randy = 1.5-3.*myrand.u();
    
    double x = (x0 + xa)*0.5 + randx;
    double y = (y0 + ya)*0.5 + randy;

    Polygon polycopy = poly;
    std::vector<Polygon> result;
    scoop_poly(polycopy, Point(x,y), result);
    
    size_t i = 0;
    for (auto &r : result)
    {
      plot_poly( r, "test_polygon2_" + std::to_string(j) + "_" + std::to_string(i) );
      ++i;
    }
  }
}


void test_polygon3()
{
  const double x0 = -0.22;
  const double y0 =  0.11;
  const auto a = 1 / sqrt(2);
  const double xa = x0+a;
  const double ya = y0+a;
  const Polygon poly =
  {
    Segment(Point(x0,y0), Point(xa,y0), XPLUS,  Point(x0,y0)),
    Segment(Point(xa,y0), Point(xa,ya), YPLUS,  Point(xa,y0)),
    Segment(Point(xa,ya), Point(x0,ya), XMINUS, Point(xa,ya)),
    Segment(Point(x0,ya), Point(x0,y0), YMINUS, Point(x0,ya))
  };
  
  const double maxbite = 0.01;
  const double minbite = 0.000001;
  // const double dbite = 0.001;
  const double dbite = 0.004;

  // const double dstep = 0.098;
  const double dstep = 0.298;

  // bottom
  {
    size_t j = 0;
    for (double y=y0-1.+minbite; y<=y0-1.+maxbite; y+=dbite)
    {
      for (double x=x0; x<xa; x+=dstep)
      {
        Polygon polycopy = poly;
        std::vector<Polygon> result;
        scoop_poly(polycopy, Point(x,y), result);
        
        size_t i = 0;
        for (auto &r : result)
        {
          plot_poly( r, "test_polygon3_bottom_" + std::to_string(j) + "_" + std::to_string(i) );
          ++i;
        }
        ++j;
      }
    }
  }
  
  // top
  {
    size_t j = 0;
    for (double y=ya+1.-minbite; y>=ya+1.-maxbite; y-=dbite)
    {
      for (double x=x0; x<xa; x+=dstep)
      {
        Polygon polycopy = poly;
        std::vector<Polygon> result;
        scoop_poly(polycopy, Point(x,y), result);
        
        size_t i = 0;
        for (auto &r : result)
        {
          plot_poly( r, "test_polygon3_top_" + std::to_string(j) + "_" + std::to_string(i) );
          ++i;
        }
        ++j;
      }
    }
  }

  // left
  {
    size_t j = 0;
    for (double x=x0-1.+minbite; x>=x0-1.+maxbite; x+=dbite)
    {
      for (double y=y0; y<ya; y+=dstep)
      {
        Polygon polycopy = poly;
        std::vector<Polygon> result;
        scoop_poly(polycopy,  Point(x,y), result);
        
        size_t i = 0;
        for (auto &r : result)
        {
          plot_poly( r, "test_polygon3_left_" + std::to_string(j) + "_" + std::to_string(i) );
          ++i;
        }
        ++j;
      }
    }
  }
  
  // right
  {
    size_t j = 0;
    for (double x=xa+1.-minbite; x<=xa+1.-maxbite; x-=dbite)
    {
      for (double y=y0; y<ya; y+=dstep)
      {
        Polygon polycopy = poly;
        std::vector<Polygon> result;
        scoop_poly(polycopy, Point(x,y), result);
        
        size_t i = 0;
        for (auto &r : result)
        {
          plot_poly( r, "test_polygon3_right_" + std::to_string(j) + "_" + std::to_string(i) );
          ++i;
        }
        ++j;
      }
    }
  }
  
  
  // right
  {
    size_t j = 0;
    for (double y=ya+1.-minbite; y<=ya+1.-maxbite; y-=dbite)
    {
      for (double x=x0; x<xa; x+=dstep)
      {
        Polygon polycopy = poly;
        std::vector<Polygon> result;
        scoop_poly(polycopy, Point(x,y), result);
        
        size_t i = 0;
        for (auto &r : result)
        {
          plot_poly( r, "test_polygon3_top_" + std::to_string(j) + "_" + std::to_string(i) );
          ++i;
        }
        ++j;
      }
    }
  }

}

bool disk_accept(Point disk, std::vector<Point> &bites)
{
  for (auto &circ : bites)
  {
    if (!circle_outside_circle(disk, circ))
    {
      return false;
    }
  }
  return true;
}

bool disk_covers_polygon(Point disk,
                         const Polygon &poly)
{
  for (auto &s : poly)
  {
    if (point_outside_circle(s.p0,disk))
      return false;
  }
  return true;
}

bool disk_covers_a_polygon(Point disk,
                           const std::vector<Polygon> &polys)
{
  for (auto &poly : polys)
  {
    if (disk_covers_polygon(disk, poly))
    {
      return true;
    }
  }
  return false;
}

// generate the next disk, give the prior disks and the remaining void polygons
bool next_bite(double max_darts,
               std::vector<Point> &bites,
               const std::vector<Polygon> &polys,
               MyRand &myrand,
               double x0, double xa, double y0, double ya,
               bool try_partial,
               bool try_outside_box)
{
  // we only allow circles outside all prior circles, the scoop code asserts otherwise
  for (size_t darts=0; darts<max_darts; ++darts)
  {
    // random numbers
    const auto randx = (2*myrand.u() - 1.); // [-1,1]
    const auto randy = (2*myrand.u() - 1.); // [-1,1]
    
    // from center
    const double x_squarecenter = (x0+xa)*0.5;
    const double y_squarecenter = (y0+ya)*0.5;
    const double maxd = (1.+(xa-x0)*0.5); // max x or y distance from center of square to center of circle

    const auto x = x_squarecenter + maxd * randx;
    const auto y = y_squarecenter + maxd * randy;
    
    // from sides, exclude points inside the square itself, uninteresting test
    // flawed, excludes any x *or* any y coinciding with the box
    //    const auto x = (randx>0 ? xa : x0 ) + randx;
    //    const auto y = (randy>0 ? ya : y0 ) + randy;
    
    Point disk( Point(x,y) );
    
    if (try_outside_box && is_point_in_box( disk, x0, xa, y0, ya))
    {
      continue;
    }
    
    if (!disk_accept(disk, bites))
    {
      continue;
    }
    
    if (try_partial && disk_covers_a_polygon(disk, polys))
    {
      continue;
    }
    // it's a good bite, take it
    bites.emplace_back( disk );
    return true;
  }
  // couldn't find a disk
  return false;
}

bool next_bite(std::vector<Point> &bites,
               const std::vector<Polygon> &polys,
               MyRand &myrand,
               double x0, double xa, double y0, double ya)
{
  // 4000, 2000, 10000
  if (next_bite(4000, bites, polys, myrand, x0, xa, y0, ya, true, true))
    return true;
  if (next_bite(1000, bites, polys, myrand, x0, xa, y0, ya, false, true))
    return true;
  if (next_bite(1000, bites, polys, myrand, x0, xa, y0, ya, false, false))
    return true;
  return false;
}

void default_square_poly(double &x0, double &y0, double &xa, double &ya, double &a, Polygon &square)
{
  x0 = 0;
  y0 = 0;
  a = 1 / sqrt(2);
  xa = x0+a;
  ya = y0+a;
  square =
  {
    Segment(Point(x0,y0), Point(xa,y0), XPLUS,  Point(x0,y0)),
    Segment(Point(xa,y0), Point(xa,ya), YPLUS,  Point(xa,y0)),
    Segment(Point(xa,ya), Point(x0,ya), XMINUS, Point(xa,ya)),
    Segment(Point(x0,ya), Point(x0,y0), YMINUS, Point(x0,ya))
  };
}

void random_bites(std::vector<Polygon> &polys, std::vector<Point> &bites, int max_bites, MyRand &myrand)
{
  double x0, y0, xa, ya, a;
  Polygon square;
  default_square_poly(x0, y0, xa, ya, a, square );
  
  polys = {square};
  
  std::vector<Polygon> result;
  for (size_t b = 0; b < max_bites; ++b)
  {
    if (!next_bite(4000, bites, polys, myrand, x0, xa, y0, ya, true, true))
    {
      break;
    }
    auto &bite = bites.back();
    
    for (auto &poly : polys)
    {
      scoop_poly(poly, bite, result);
    }
    
    // quit if the square is covered
    if (polys.empty())
    {
      bites.pop_back();
      break;
    }
    
    // keep going if there's something
    std::swap(result,polys);
    result.clear();
  }
}

// iteratively cut
void test_polygon4()
{
  // parameters
  const size_t num_tests = 1000000; // 100;
  const size_t num_bites = 15; // E-MPS paper lemma 3 proved at most 15 can touch the inner square
  // E-MPS paper lemma 9 proved <= 9 can actually contribute to a non-empty void polygon, and 8 is probably tight.
  
  bool do_plots = false;
    
  MyRand myrand(1334324123);
  
  double x0, y0, xa, ya, a;
  Polygon square;
  default_square_poly(x0, y0, xa, ya, a, square );

  for (size_t t=0; t<num_tests; ++t)
  {
    std::vector<Polygon> polys = {square};
    std::vector<Point> bites;
    
    std::vector<Polygon> result;
    for (size_t b = 0; b < num_bites; ++b)
    {
      if (!next_bite(bites, polys, myrand, x0, xa, y0, ya))
        break;
      auto &bite = bites.back();
      
      for (auto &poly : polys)
      {
        scoop_poly(poly, bite, result);
      }
      std::swap(result,polys);
      result.clear();
      
      // plot collection of bites so far
      if (do_plots)
      {
        std::string fname = "test_polygon4_test" + std::to_string(t) + "_bite" + std::to_string(b) + "_";        
        plot_square_biters( Point(x0,y0), Point(xa,ya), bites, b, fname );
        
        // plot each polygon result, in square
        size_t i = 0;
        for (auto &poly : polys)
        {
          plot_poly(poly, fname + "poly" + std::to_string(i),
               Point(x0,y0), Point(xa,ya) );
          ++i;
        }
      }
      
      // quit if the square is covered
      if (polys.empty())
      {
        break;
      }
      
    }
  }
}

void test_scoop_circle_circle0()
{
  const double x0 = 0.;
  const double y0 = 0.;
  const auto a = 1 / sqrt(2);
  const double xa = x0+a;
  const double ya = y0+a;
  const Polygon poly =
  {
    Segment(Point(x0,y0), Point(xa,y0), XPLUS,  Point(x0,y0)),
    Segment(Point(xa,y0), Point(xa,ya), YPLUS,  Point(xa,y0)),
    Segment(Point(xa,ya), Point(x0,ya), XMINUS, Point(xa,ya)),
    Segment(Point(x0,ya), Point(x0,y0), YMINUS, Point(x0,ya))
  };
  
  Polygon polycopy = poly;
  std::vector<Polygon> result0;
  Point d0( Point(x0+a*0.5, y0 -1. + 0.55*a) );
  scoop_poly(polycopy, d0, result0);
  
  Point d1( Point(x0+a*0.5, ya +1. - 0.55*a) );
  size_t j = 0;
  for (auto &r0 : result0)
  {
    std::vector<Polygon> result1;
    scoop_poly(r0, d1, result1);
    
    size_t i = 0;
    for (auto &r1 : result1)
    {
      plot_poly( r1, "test_scoop_circle_circle0_" + std::to_string(j) + "_" + std::to_string(i) );
      ++i;
    }
    ++j;
  }
}


void test_scoop_circle_circle1()
{
  const double x0 = 0.;
  const double y0 = 0.;
  const auto a = 1 / sqrt(2);
  const double xa = x0+a;
  const double ya = y0+a;
  const Polygon poly =
  {
    Segment(Point(x0,y0), Point(xa,y0), XPLUS,  Point(x0,y0)),
    Segment(Point(xa,y0), Point(xa,ya), YPLUS,  Point(xa,y0)),
    Segment(Point(xa,ya), Point(x0,ya), XMINUS, Point(xa,ya)),
    Segment(Point(x0,ya), Point(x0,y0), YMINUS, Point(x0,ya))
  };
  
  Polygon polycopy = poly;
  std::vector<Polygon> result0;
  // works
  // Circle d1( Point(x0, y0 -1. + 0.55*a) );
  // Circle d0( Point(x0, ya +1. - 0.55*a) );
  // fails
  Point d0( Point(x0, y0 -1. + 0.55*a) );
  Point d1( Point(x0, ya +1. - 0.55*a) );
  scoop_poly(polycopy, d0, result0);
  
  size_t j = 0;
  for (auto &r0 : result0)
  {
    std::vector<Polygon> result1;
    plot_poly( r0, "test_scoop_circle_circle1_" + std::to_string(j) + "_before", Point(x0,y0), Point(xa,ya) );
    scoop_poly(r0, d1, result1);
    
    size_t i = 0;
    for (auto &r1 : result1)
    {
      plot_poly( r1, "test_scoop_circle_circle1_" + std::to_string(j) + "_" + std::to_string(i) );
      ++i;
    }
    ++j;
  }
}

void test_scoop_circle_circle2()
{
  const double x0 = 0.;
  const double y0 = 0.;
  const auto a = 1 / sqrt(2);
  const double xa = x0+a;
  const double ya = y0+a;
  const Polygon poly =
  {
    Segment(Point(x0,y0), Point(xa,y0), XPLUS,  Point(x0,y0)),
    Segment(Point(xa,y0), Point(xa,ya), YPLUS,  Point(xa,y0)),
    Segment(Point(xa,ya), Point(x0,ya), XMINUS, Point(xa,ya)),
    Segment(Point(x0,ya), Point(x0,y0), YMINUS, Point(x0,ya))
  };
  
  Polygon polycopy = poly;
  std::vector<Polygon> result0;
  Point d0( Point(xa, y0 -1. + 0.55*a) );
  scoop_poly(polycopy, d0, result0);
  
  Point d1( Point(xa, ya +1. - 0.55*a) );
  size_t j = 0;
  for (auto &r0 : result0)
  {
    std::vector<Polygon> result1;
    scoop_poly(r0, d1, result1);
    
    size_t i = 0;
    for (auto &r1 : result1)
    {
      plot_poly( r1, "test_scoop_circle_circle2_" + std::to_string(j) + "_" + std::to_string(i) );
      ++i;
    }
    ++j;
  }
}

void test_scoop_degeneracy_1()
{
  const auto a = sqrt(0.5); // 1 / sqrt(2);
  const auto a2 = 0.5*a;
  const double x0 = a;
  const double y0 = a;
  const double xa = x0+a;
  const double ya = y0+a;
  const Polygon poly =
  {
    Segment(Point(x0,y0), Point(xa,y0), XPLUS,  Point(x0,y0)),
    Segment(Point(xa,y0), Point(xa,ya), YPLUS,  Point(xa,y0)),
    Segment(Point(xa,ya), Point(x0,ya), XMINUS, Point(xa,ya)),
    Segment(Point(x0,ya), Point(x0,y0), YMINUS, Point(x0,ya))
  };

  Point d0( Point(a2, a2) );
  Point d1( Point(a2, a2+a*2.) );
  Point d2( Point(a2+a*2., a2) );
  Point d3( Point(a2+a*2., a2+a*2.) );
  
  std::vector<Polygon> poly0, poly1, poly2, poly3;
  scoop_poly(poly,  d0, poly0);
  scoop_poly(poly0.front(), d1, poly1);
  scoop_poly(poly1.front(), d2, poly2);
  scoop_poly(poly2.front(), d3, poly3, true); // this step fails to produce a good poly

  plot_poly( poly, "test_scoop_degeneracy_1-poly" );
  plot_poly( poly0.front(), "test_scoop_degeneracy_1-poly0" );
  plot_poly( poly1.front(), "test_scoop_degeneracy_1-poly1" );
  plot_poly( poly2.front(), "test_scoop_degeneracy_1-poly2" );
  plot_poly( poly3.front(), "test_scoop_degeneracy_1-poly3" );
}
