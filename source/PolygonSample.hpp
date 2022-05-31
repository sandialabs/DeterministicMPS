//
//  PolygonSample.hpp
//  AvoidReject
//
//  Created by Mitchell, Scott A on 10/21/21.
//  Copyright © 2021 Mitchell, Scott A. All rights reserved.
//

#ifndef PolygonSample_hpp
#define PolygonSample_hpp

#include "MyRand.hpp"
#include "Geometry.hpp"
#include "Triangulate.hpp"

// triangle abc area
double triangle_area(const Point &a, const Point &b, const Point &c);

// return random point P in triangle abc, given uniform random u and v in [0,1]
void sample_triangle(const Point &a, const Point &b, const Point &c, double u, double v, Point &P);

double polygon_area(const Chocks &chocks, const Loop &loop, const Triangles &triangles);

// return index of chock or triangle, uniform random by area given sub_area = u * polygon_area, u uniform random [0,1]
// -1 if not picked
void pick_subregion(const Chocks &chocks, const Loop &loop, const Triangles &triangles, double sub_area,
                    int &ci, int &ti);

// return random point p from polygon
//   prior to calling, shave chocks and triangulate loop
// poly_area can be provided as a speedup, otherwise we'll calculate it from scratch
void sample_polygon(const Chocks &chocks, const Loop &loop, const Triangles &triangles,
                    Point &p,
                    MyRand &myrand, double poly_area = -1.);

#endif /* PolygonSample_hpp */
