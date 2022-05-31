//
//  TestTriangulate.hpp
//  AvoidReject
//
//  Created by Mitchell, Scott A on 10/20/21.
//  Copyright © 2021 Mitchell, Scott A. All rights reserved.
//

#ifndef TestTriangulate_hpp
#define TestTriangulate_hpp

#include "Geometry.hpp"
#include "Triangulate.hpp"

class MyRand;

// void test_shave(const Polygon &poly, Chocks &chocks, Loop &loop, Point bbox_lo, Point bbox_hi, std::string fname);

void test_shave_0(); // two little voids
void test_shave_1(); // one big void
void test_shave_2(); // 100 random voids, of varying number of sides
void test_shave_3(); // 1,000,000 random voids, of varying number of sides, no output

// triangulate same voids as shave
void test_triangulate_0();
void test_triangulate_1();
void test_triangulate_2();
void test_triangulate_3();

// unit tests of functions used by triangulations
void test_angle_measure();
void test_circum();


// stress test, moving two disks close together until they overlap, ensure the triangulation always succeeds
// vertical, horizontal, and at odd angles, halfing the distance between each iteration
void test_grazing_0(); // id 5 6
// same idea but moving close to the square sides
void test_grazing_1(); // id 7 8
void test_grazing_2(); // id 9 10

// return the polys that result from up to max_bites from a square
void random_bites(std::vector<Polygon> &polys, std::vector<Point> &bites, int max_bites, MyRand &myrand);

// figures for paper
void figure_trim();
void figure_triangles(); 


#endif /* TestTriangulate_hpp */
