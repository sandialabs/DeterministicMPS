//
//  PolygonSample.cpp
//  AvoidReject
//
//  Created by Mitchell, Scott A on 10/21/21.
//  Copyright © 2021 Mitchell, Scott A. All rights reserved.
//

#include "PolygonSample.hpp"

#include <iostream>

#include "ChockSample.hpp"

void sample_triangle(const Point &a, const Point &b, const Point &c, double u, double v, Point &P)
{
  if (u + v > 1.)
  {
    u = 1. - u;
    v = 1. - v;
  }
  
  // A + AB * p + BC * q
  P.x = a.x + (b.x - a.x) * u + (c.x - a.x) * v;
  P.y = a.y + (b.y - a.y) * u + (c.y - a.y) * v;
}

double polygon_area(const Chocks &chocks, const Loop &loop, const Triangles &triangles)
{
  double area=0.;
  for (auto &chock : chocks)
  {
    area += chock.area();
  }
  for (size_t t = 0; t < triangles.size(); t+=3)
  {
    area += triangle_area( loop[triangles[t]], loop[triangles[t+1]], loop[triangles[t+2]]);
  }
  return area;
}

void pick_subregion(const Chocks &chocks, const Loop &loop, const Triangles &triangles, double sub_area,
                    int &ci, int &ti)
{
  if (chocks.empty() && triangles.empty())
  {
    std::cout << "ERROR: pick_subregion failed. No chocks or triangles to choose from." << std::endl;
    ci=-1;
    ti=-1;
    return;
  }
  double area=0.;
  for (size_t i = 0; i < chocks.size(); ++i)
  {
    auto &chock = chocks[i];
    area += chock.area();
    if (area >= sub_area)
    {
      ci = (int) i;
      ti = -1;
      return;
    }
  }
  for (size_t t = 0; t < triangles.size(); t+=3)
  {
    area += triangle_area( loop[triangles[t]], loop[triangles[t+1]], loop[triangles[t+2]]);
    if (area >= sub_area)
    {
      ti = (int) t;
      ci = -1;
      return;
    }
  }
  // error, sub_area was too large
  std::cout << "ERROR: pick_subregion failed. Total area is " << area <<
  " but larger sub_area " << sub_area << " requested." << std::endl;
  // one could use the last triangle or chock anyway
  if (triangles.empty())
  {
    ci = (int) chocks.size()-1;
    ti = -1;
  }
  else
  {
    ti = (int) triangles.size() - 1;
    ci = -1;
  }
}

void sample_polygon(const Chocks &chocks, const Loop &loop, const Triangles &triangles,
                    Point &p,
                    MyRand &myrand, double poly_area )
{
  if (poly_area<0.)
  {
    poly_area = polygon_area(chocks, loop, triangles);
  }
  
  const double region_u = myrand.u();
  const double sub_area = poly_area * region_u;
  
  int ci, ti;
  pick_subregion(chocks, loop, triangles, sub_area, ci, ti);
  
  // random point in subregion
  const double u = myrand.u();
  const double v = myrand.u();
  
  // point from chock
  if (ci >= 0)
  {
    const Chock &chock = chocks[ci];
    
    double phi_s, r_s;
    sample_chock( u, v, chock.phi(), phi_s, r_s);
    chock_to_cartesian(chock.c, chock.t, chock.q, phi_s, r_s, p.x, p.y);
  }
  // point from triangle
  else if (ti >= 0)
  {
    sample_triangle(loop[triangles[ti]], loop[triangles[ti+1]], loop[triangles[ti+2]], u, v, p);
  }
  else
  {
    if (!loop.empty())
    {
      p = loop.front();
    }
    else if (!chocks.empty())
    {
      p = chocks[0].q;
    }
    else
    {
      std::cout << "ERROR: sample_polygon failed because polygon is empty." << std::endl;
      assert(0);
    }
  }
}

void sample_polygons(const std::vector< Chocks > &vchocks,
                     const std::vector<Loop> &vloop,
                     const std::vector< Triangles > &triangles,
                     std::vector<double> vpoly_area,
                     Point &p,
                     MyRand &myrand);
