//
//  Triangulate.hpp
//  AvoidReject
//
//  Created by Mitchell, Scott A on 10/19/21.
//  Copyright © 2021 Mitchell, Scott A. All rights reserved.
//

#ifndef Triangulate_hpp
#define Triangulate_hpp

#include "Geometry.hpp"

// Polygons have segments, which are defined by circles or sides of squares
// Loops are composed only of straight segments, defined by a chain of vertices, in CCW order, not repeated.
// Chocks are defined by three points: c, t, q; circle center, circle tangent, and apex.
// Triangles are indices into the loop of which vertices define triangles.

typedef std::vector<Point> Loop;
typedef std::vector<int> Triangles;

void plot_loop(const Loop &loop, std::ostream &fout, bool do_labels = true, bool do_fill = false );
void print_loop(const Loop &loop);

void seg_circle_to_chock(const Polygon &poly, size_t si, Chocks &chocks, Loop &loop);
void seg_square_to_loop(const Polygon &poly, size_t si, Loop &loop);

// decompose a polygon into chocks and one loop
void shave_chocks(const Polygon &poly, Chocks &chocks, Loop &loop);

// remove repeated points
void simplify_loop(Loop &loop);

// decompose a loop into triangles
// return error code: 0 = no error, 1 = warning, 2 = error
// makes prettier pictures
int triangulate_loop( const Loop &loop, Triangles &triangles); // with care taken for triangle shape
// more robust for nearly colinear or coincident points
int triangulate_loop_simple( const Loop &loop, Triangles &triangles); // only care about triangle orientation and empty

void plot_triangles(const Loop &loop, const Triangles &triangles, std::ostream &fout, bool do_labels = true);
void plot_triangulation(const Loop &loop, const Triangles &triangles, bool do_plot_box, Point bbox_lo, Point bbox_hi, std::string fname,
                        bool do_labels = true);

// exposed for testing only
// monotonically increasing function of angle at b
double measure_angle(const Point &a, const Point &b, const Point &c);
void circumcenter(const Point &a, const Point &b, const Point &c, Point &o, double &R );
bool is_empty_circle_b(const Point &a, const Point &b, const Point &c, const Point &m,
                       const Point &o, double R);

#endif /* Triangulate_hpp */
