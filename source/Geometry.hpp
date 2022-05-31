//
//  Geometry.hpp
//  AvoidReject
//
//  Created by Mitchell, Scott A on 7/28/21.
//  Copyright © 2021 Mitchell, Scott A. All rights reserved.
//

#ifndef Geometry_hpp
#define Geometry_hpp

#include <cmath>
#include <vector>
#include <string>
#include <ostream>

// =================== interface
// ==== point
class Point // or vector
{
public:
  double x, y;
  
  Point(double x0, double y0) : x(x0), y(y0) {}
  Point() : x(0.), y(0.) {}
  Point( const Point &other ) : x(other.x), y(other.y) {}

  Point& operator=(const Point& other);
  Point& operator+=(const Point& rhs);
  friend Point operator+(Point lhs, const Point& rhs);
  Point& operator-=(const Point& rhs);
  friend Point operator-(Point lhs, const Point& rhs);
  Point& operator*=(double m);
  Point  operator*(double m) const;
  Point  operator- ();

  bool   operator==(const Point &q) const;
  bool   operator!=(const Point &q) const;
};

double cross_product( const Point &p, const Point &q);
double dot_product( const Point &p, const Point &q);
double norm( const Point &p );
double normsquared( const Point &p );
double atan2( const Point &p );
// Return true if the shorter rotation of p onto q is clockwise
// p and q are vectors, assume shared starting point.
bool clockwise( const Point &p, const Point &q);
bool clockwise_strict( const Point &p, const Point &q); // not colinear

bool is_point_in_box( const Point &p, double x0, double xa, double y0, double ya );

void print_point(const Point &p);

// ==== circle

// used in early versions for testing chock sampling with variable radii
//
//class CircleWithRadius
//{
//public:
//  Point c;
//  const double r=1.;
//
//  CircleWithRadius() : c(), r(1.) {}
//  CircleWithRadius(const Point &p) : c(p), r(1.) {}
//  CircleWithRadius(double rr) : c(), r(rr) {}
//
//  // bool outside(Point p) const;
//
//  // we intersect with axis-aligned segments, and other Circles only
//};

class Circle
{
public:
  Point c;
  // static const double r=1.;
  
  Circle() : c() {}
  Circle(const Point &p) : c(p) {}
  
  bool outside(Point p) const;
  // Symmetric. true if circ center is outside this circle, and vice versa. Symmetric.
  bool outside(const Circle &circ) const;
  // a sample is OK if it is not_inside any disk
  bool not_inside(Point p) const;

  // we intersect with axis-aligned segments, and other Circles only
  
  void plot_circle( std::ostream &fout ) const;
};

// true if p is outside the unit circle with center c
bool point_outside_circle(const Point &p, const Point &c);

// Symmetric. true if circ center c is outside d, and vice versa. Symmetric.
bool circle_outside_circle(const Point &c, const Point &d);

// true if p is not inside the unit circle at center c
bool point_not_inside_circle(const Point &p, const Point &c);

void plot_circle( const Point &c, std::ostream &fout );

// intersect two circles with centers d0 and d1, assuming centers are outside the other circle
void intersect_circle(const Point &d0, const Point &d1, bool &disjoint, Point &n0, Point &n1);

// intersect circle with vertical line at x=x
void intersect_x(const Circle &d, double x, bool &disjoint, Point &n0, Point &n1);

// intersect circle with horizontal line at y=y
void intersect_y(Circle d, double y, bool &disjoint, Point &n0, Point &n1);

// true if p to q in the short direction is clockwise about the circle center d
bool clockwise(const Point &d, const Point &p, const Point &q);

// return the angles of p and q from the circle center d, with pangle < 2pi and qangle > pangle
void angles(const Point &d, const Point &p, const Point &q, double &pangle, double &qangle);

// rotate angle a by 2pi to try to make it lie between b and c, if possible
// i.e. change a to be in [b,c], by adding multiples of 2pi to a
// return true if success; else a is still outside that interval
bool make_between( double &a, double b, double c);

// convert point at angle theta on circle circ to cartesian x,y
void theta_to_cartesian( const Point &circ, double theta, Point &p);
void theta_to_cartesian( const Circle &circ, double theta, double r, Point &p);

// ==== chock
struct Chock
{
  // c = center of circle, assumed radius 1
  // t = tangent point on circle
  // q = point outside circle at apex of chock
  Point c, t, q;
  
  Chock( const Point &cc, const Point &tt, const Point &qq ) : c(cc), t(tt), q(qq) {}
  
  double phi() const;
  double area() const;
};
typedef std::vector<Chock> Chocks;

// convert a point inside the chock given as (phi, r) relative to c and t, into plane coordinates (x,y)
void chock_to_cartesian( const Point&c, const Point&t, const Point&q, double phi, double r, double &x, double &y );

// phi as a positive angle < pi/2
double chock_phi( const Point&c, const Point&t, const Point&q );

// area of chock for given phi
double chock_area( double phi );

// extent of radius for given phi
double chock_r( double phi);

void plot_chock(const Point &c, const Point &t, const Point &q, std::ostream &fout, const bool do_center = true);
void plot_chocks(const Chocks &chocks, std::ostream &fout, const bool do_center = true);
void plot_full_circles(const Chocks &chocks, std::ostream &fout);

// ==== segments of polygons
enum CurveType {NADA, CIRCLE, XPLUS, XMINUS, YPLUS, YMINUS};
// xplus means horizontal and x0 < x1 is the ccw orientation

// not making Segment a virtual class so we don't have to deal with pointers and memory allocation
class Segment
{
public:
  Point p0, p1;
  // Curve type is derived: horizontal line, vertical line, circle

  Point dp; // circle if this is a CIRCLE, else x or y intercept for lines
  CurveType ctype = NADA; // NADA indicates null segment

  Segment() {}
  Segment(Point q0, Point q1) : p0(q0), p1(q1) {}
  Segment(Point q0, Point q1, CurveType dtype, Point d) : p0(q0), p1(q1), ctype(dtype), dp(d) {}

  void plot_segment( std::ostream &fout ) const;

  // true if this segment is a segment of the circle centered at c
  bool is_circle(const Point &d) const;
  
  bool   operator==(const Segment &q) const;
};

// ==== polygons

typedef std::vector<Segment> Polygon; // segments in ccw order, segs[0].s1 equivalent to segs[1].s0, etc.

// return number of segments added to poly
//   set inside true if seg was entirely inside disk d
//   set outside true if seg was entirely outside disk d
int scoop_segment(const Segment &seg, const Point &d, Polygon &poly, bool &inside, bool &outside);

// construct polygon for grid square with lower corner (x0,y0)
Polygon square_to_poly(double x0, double y0);

// return bounding box of points of polygon
void bounding_box(const Polygon &poly, Point &bbox_lo, Point &bbox_hi);

typedef std::vector<Polygon> ScoopWorkspace;

// scoop circle d out of polygon poly. Result is in some set of polygons
void scoop_poly(const Polygon &poly, const Point &d, std::vector<Polygon> &result, bool do_debug=false);
// scoop circle d out of all polys, replace polys with the resulting polygons
//   workspace is provided for efficiency, swapping
void scoop_polys(std::vector<Polygon> &polys, const Point &d, std::vector<Polygon> &workspace, bool check_covered = false);

// true if the poly is completely covered by the disk d
bool poly_covered(const Polygon &poly, const Point &d);


void plot_poly( const Polygon &poly, std::string fname );
void plot_poly( const Polygon &poly, std::string fname, Point bbox_lo, Point bbox_hi );
void print_poly(const Polygon &poly);

// plot a square and a bunch of circles, nothing clipped
void plot_square_biters(Point square_lo, Point square_hi, const std::vector<Point> &bites, size_t b, std::string fname);

// set scale, linewidth, etc.
void plot_preamble( std::ofstream &fout );

void plot_bounding_box(Point bbox_lo, Point bbox_hi, bool expand_box, std::ofstream &fout);



// ============ implementations

inline
double cross_product( const Point &p, const Point &q)
{
  return p.x * q.y - q.x * p.y;
}

inline
double dot_product( const Point &p, const Point &q)
{
  return p.x * q.x + q.y * p.y;
}

inline
double norm( const Point &p )
{
  return sqrt(p.x * p.x + p.y * p.y);
}

inline
double normsquared( const Point &p )
{
  return p.x * p.x + p.y * p.y;
}

inline
double atan2( const Point &p )
{
  return atan2(p.y,p.x);
}

inline
bool clockwise( const Point &p, const Point &q)
{
  return (cross_product(p,q)<=0.);
}
inline
bool clockwise_strict( const Point &p, const Point &q)
{
  return (cross_product(p,q)<0.);
}

inline
Point& Point::operator=(const Point& other)
{
  x=other.x;
  y=other.y;
  return *this;
}

inline
Point& Point::operator+=(const Point& rhs)
{
  x+=rhs.x;
  y+=rhs.y;
  return *this;
}

inline
Point operator+(Point lhs, const Point& rhs)
{
  lhs += rhs;
  return lhs;
}

inline
Point& Point::operator-=(const Point& rhs)
{
  x-=rhs.x;
  y-=rhs.y;
  return *this;
}

inline
Point operator-(Point lhs, const Point& rhs)
{
  lhs -= rhs;
  return lhs;
}

inline
Point& Point::operator*=(double m)
{
  x*=m;
  y*=m;
  return *this;
}

inline
Point Point::operator*(double m) const
{
  Point p( *this );
  p *= m;
  return p;
}

inline
bool Point::operator==(const Point &q) const
{
  return ((x == q.x) && (y == q.y));
}

inline
bool Point::operator!=(const Point &q) const
{
  return ((x != q.x) || (y != q.y));
}

inline
Point Point::operator- ()
{
  return Point(-x,-y);
}


//== circle

inline
bool outside(const Circle &d, Point p)
{
  return (normsquared( p-d.c ) >= 1. /* d.r*d.r */ );
}

inline
void angles(const Point &d, const Point &p, const Point &q, double &pangle, double &qangle)
{
  pangle = atan2( p - d );
  qangle = atan2( q - d );
  if (qangle < pangle)
    qangle += 2.*M_PI;
}
          
//== chock
inline
double chock_area(double phi)
{
  return (tan(phi) - phi)/2.;
}

inline
double chock_r(double phi)
{
  return 1./cos(phi); // secant(phi)
}

inline
bool clockwise(const Point &d, const Point &p, const Point &q)
{
  return clockwise(p-d,q-d);
}

inline
bool Segment::is_circle(const Point &d) const
{
  return ctype == CIRCLE && dp == d;
}

inline
bool Circle::outside(Point p) const
{
  return (normsquared( p - c ) > 1.);
}
inline
bool Circle::not_inside(Point p) const
{
  return (normsquared( p - c ) >= 1.);
}

inline
bool Circle::outside(const Circle &circ) const
{
  return ( outside(circ.c) && circ.outside(c) );
}


// true if p is outside the unit circle with center c
inline
bool point_outside_circle(const Point &p, const Point &c)
{
  return (normsquared(p-c)>1.);
}

// Symmetric. true if circ center c is outside d, and vice versa. Symmetric.
inline
bool circle_outside_circle(const Point &c, const Point &d)
{
  return (normsquared(c-d)>1.);
}

// true if p is not inside the unit circle at center c
inline
bool point_not_inside_circle(const Point &p, const Point &c)
{
  return (normsquared(p-c)>=1.);
}


inline
bool is_point_in_box( const Point &p, double x0, double xa, double y0, double ya )
{
  return (p.x>=x0 && p.x<=xa && p.y>=y0 && p.y<=ya);
}

// convert point at angle theta on circle circ to cartesian x,y
inline
void theta_to_cartesian( const Point &circ, double theta, Point &p)
{
  p.x = cos(theta) + circ.x;
  p.y = sin(theta) + circ.y;
}

inline
void theta_to_cartesian( const Circle &circ, double theta, double r, Point &p)
{
  p.x = r*cos(theta) + circ.c.x;
  p.y = r*sin(theta) + circ.c.y;
}

inline
double chock_phi( const Point&c, const Point&t, const Point&q )
{
  return atan(norm(q-t));
}

inline
double Chock::phi() const
{
  return chock_phi(c,t,q);
}

inline
double Chock::area() const
{
  return chock_area(phi());
}

inline
bool Segment::operator==(const Segment &q) const
{
  return
  ctype == q.ctype &&
  p0 == q.p0 &&
  p1 == q.p1 &&
  dp == q.dp;
}

#endif /* Geometry_hpp */
