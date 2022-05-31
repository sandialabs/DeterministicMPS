//
//  TestGeometry.hpp
//  AvoidReject
//
//  Created by Mitchell, Scott A on 8/20/21.
//  Copyright © 2021 Mitchell, Scott A. All rights reserved.
//

#ifndef TestGeometry_hpp
#define TestGeometry_hpp

#include "Geometry.hpp"
#include "MyRand.hpp"

void test_intersect_circle();
void test_intersect_x();
void test_intersect_y();

// test scooping a square with one disk
void test_polygon0(); // 1 of each of known cases
void test_polygon1(); // uniform spacing
void test_polygon2(); // random placement
void test_polygon3(); // barely-bit
void test_polygon4(); // iterative bites

void test_scoop_circle_circle0();
void test_scoop_circle_circle1();
void test_scoop_circle_circle2();
#endif /* TestGeometry_hpp */

void test_scoop_degeneracy_1();

void default_square_poly(double &x0, double &y0, double &xa, double &ya, double &a, Polygon &square);
void random_bites(std::vector<Polygon> &polys, std::vector<Point> &bites, int max_bites, MyRand &myrand);

