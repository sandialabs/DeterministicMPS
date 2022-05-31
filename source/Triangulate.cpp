//
//  Triangulate.cpp
//  AvoidReject
//
//  Created by Mitchell, Scott A on 10/19/21.
//  Copyright © 2021 Mitchell, Scott A. All rights reserved.
//
// For triangulating polygons into chocks and triangles
//
#include "Triangulate.hpp"

#include <iostream>
#include <fstream>
#include <numeric>
#include <algorithm>
#include <cassert>

// To ensure the loop is a simple polygon for the triangulation,
//   ensure q is strictly closer to this_center than any other center, and strictly inside box
void pull_back_q( const Polygon &poly, size_t si, Point &q )
{
  const Segment &this_seg = poly[si];
  
  for (size_t i = 0; i <poly.size(); ++i)
  {
    if (i==si)
    {
      continue;
    }
    const Segment &seg = poly[i];
    switch (seg.ctype)
    {
      case CIRCLE:
      {
        double f = 1.-__DBL_EPSILON__;
        while (normsquared(q- seg.dp) <= normsquared(q- this_seg.dp))
        {
          q = (q-this_seg.dp)*f + this_seg.dp;
          f -= __DBL_EPSILON__;
        }
      }
        break;
      case XPLUS:
      {
        const auto y = seg.dp.y;
        if (this_seg.dp.y>y)
        {
          double f = __DBL_EPSILON__;
          while (q.y <= y )
          {
            q.y = y + f;
            f += __DBL_EPSILON__;
          }
        }
      }
        break;
      case XMINUS:
      {
        const auto y = seg.dp.y;
        if (this_seg.dp.y<y)
        {
          double f = __DBL_EPSILON__;
          while( q.y >= y)
          {
            q.y = y - f;
            f += __DBL_EPSILON__;
          }
        }
      }
        break;
      case YPLUS:
      {
        const auto x = seg.dp.x;
        if (this_seg.dp.x<x)
        {
          double f = __DBL_EPSILON__;
          while (q.x >= x)
          {
            q.x = x - f;
            f += __DBL_EPSILON__;
          }
        }
      }
        break;
      case YMINUS:
      {
        const auto x = seg.dp.x;
        if (this_seg.dp.x>x)
        {
          double f = __DBL_EPSILON__;
          while(q.x <= x)
          {
            q.x = x + f;
            f += __DBL_EPSILON__;
          }
        }
      }
        break;
      default:
        break;
    }
  }
}

// To ensure the loop is a simple polygon for the triangulation,
//   ensure t is strictly closer to this_center than any other center except it's neighbor, and strictly inside box
void pull_back_t( const Polygon &poly, size_t si, Point &t )
{
  const Segment &this_seg = poly[si];
  
  for (size_t i = 0; i <poly.size(); ++i)
  {
    if (i==si || i+1==si || (si==0 && i+1== poly.size()))
    {
      continue;
    }
    const Segment &seg = poly[i];
    switch (seg.ctype)
    {
      case CIRCLE:
      {
        double f = 1.-__DBL_EPSILON__;
        while (normsquared(t- seg.dp) <= normsquared(t- this_seg.dp))
        {
          t = (t-this_seg.dp)*f + this_seg.dp;
          f -= __DBL_EPSILON__;
        }
      }
        break;
      case XPLUS:
      {
        const auto y = seg.dp.y;
        if (this_seg.dp.y>y)
        {
          double f = __DBL_EPSILON__;
          while (t.y < y )
          {
            t.y = y + f;
            f += __DBL_EPSILON__;
          }
        }
      }
        break;
      case XMINUS:
      {
        const auto y = seg.dp.y;
        if (this_seg.dp.y<y)
        {
          double f = __DBL_EPSILON__;
          while( t.y > y)
          {
            t.y = y - f;
            f += __DBL_EPSILON__;
          }
        }
      }
        break;
      case YPLUS:
      {
        const auto x = seg.dp.x;
        if (this_seg.dp.x<x)
        {
          double f = __DBL_EPSILON__;
          while (t.x > x)
          {
            t.x = x - f;
            f += __DBL_EPSILON__;
          }
        }
      }
        break;
      case YMINUS:
      {
        const auto x = seg.dp.x;
        if (this_seg.dp.x>x)
        {
          double f = __DBL_EPSILON__;
          while(t.x <= x)
          {
            t.x = x + f;
            f += __DBL_EPSILON__;
          }
        }
      }
        break;
      default:
        break;
    }
  }
}

// workspace
std::vector<double> seg_circle_to_chock_thetas;

void seg_circle_to_chock(const Polygon &poly, size_t si, Chocks &chocks, Loop &loop)
{
  const Segment &this_seg = poly[si];
  assert(this_seg.ctype==CIRCLE);
  std::vector<double> &thetas = seg_circle_to_chock_thetas;
  thetas.clear();
  thetas.reserve(poly.size());
  
  // get angles of p0 and p1
  double p0angle, p1angle;
  const auto &this_center = this_seg.dp;
  assert( clockwise( this_center, this_seg.p0, this_seg.p1) );
  angles( this_center, this_seg.p1, this_seg.p0, p1angle, p0angle);

  // optimization: fewer chocks splits if we check the max possible extent of the chock first
  // extent of chock
  assert(p0angle>=p1angle);
  auto phi = (p0angle - p1angle) * 0.5;
  const auto extended_r = chock_r(phi);
  
  const double eps = 0.001;
  const double nearly_one = 1.-eps;

  thetas.push_back(p0angle);
  
  // get angles of closest points
  for (size_t i = 0; i <poly.size(); ++i)
  {
    if (i==si)
    {
      continue;
    }
    
    // Skip trying to get this better.
    // The current test is based on drawing the circle around the chock
    // for a tighter check, instead, compute the point q of the two chocks and check
    // if q falls on the wrong side of the box line, or inside some other circle-segment's chock triangle.
    // see test_shave_2_sub55_0_triangles.ps and 81 for a test example for a box side
    // see test_shave_2_sub67_0_triangles.ps and 84, 86 for a test example with two circles
    
    bool found = false;
    double theta = 0.;
    const Segment &seg = poly[i];
    // vector x,y from circle center to closest point of exended segment.
    switch (seg.ctype)
    {
      case CIRCLE:
      {
        const double x = seg.dp.x -this_seg.dp.x;
        const double y = seg.dp.y -this_seg.dp.y;
        // no need if circles overlap
        const auto centers_distance_squared = x*x + y*y;
        if (centers_distance_squared >= 4.0 * nearly_one * nearly_one)
        {
          // circles don't overlap, or at least don't overalap a lot

          // no need if they are so far apart their chocks won't overlap
          double seg_p1angle, seg_p0angle;
          angles(seg.dp, seg.p1, seg.p0, seg_p1angle, seg_p0angle);
          auto seg_phi = (p0angle - p1angle) * 0.5;
          const auto seg_extended_r = chock_r(seg_phi);
          const auto extended_r2 = (extended_r + seg_extended_r) * (extended_r + seg_extended_r);
          if (centers_distance_squared <= extended_r2)
          {
            theta = atan2(y,x);
            found = true;
          }
        }
      }
        break;
      case XPLUS:
      {
        //if (this_seg.dp.y > seg.p0.y + 0.999)
        if ((this_seg.dp.y > seg.p0.y + nearly_one) &&   // doesn't cross the box
            (this_seg.dp.y < seg.p0.y + extended_r + eps)) // but raw chock might
        {
          theta = (M_PI * 1.5);
          found = true;
        }
      }
        break;
      case XMINUS:
      {
        if ((this_seg.dp.y < seg.p0.y - nearly_one) &&
            (this_seg.dp.y > seg.p0.y - extended_r - eps))
        {
          theta = M_PI * 0.5;
          found = true;
        }
      }
        break;
      case YPLUS:
      {
        if ((this_seg.dp.x < seg.p0.x - nearly_one) &&
            (this_seg.dp.x > seg.p0.x - extended_r - eps))
        {
          theta = 0.;
          found = true;
        }
      }
        break;
      case YMINUS:
      {
        if ((this_seg.dp.x > seg.p0.x + nearly_one) &&
            (this_seg.dp.x < seg.p0.x + extended_r + eps))
        {
          theta = M_PI;
          found = true;
        }
      }
        break;
      default:
        break;
    }
    if (found && make_between(theta,p1angle,p0angle))
    {
      thetas.push_back(theta);
    }
  }
  thetas.push_back(p1angle);
  
  // descending order
  std::sort(thetas.begin(), thetas.end(), std::greater<double>());
  
  // convert thetas to points, and add in chock-q points
  Point t0, t1, q;
  for (auto i = 0; i+1 < thetas.size(); ++i)
  {
    double theta0 = thetas[i];
    double theta1 = thetas[i+1];
    theta_to_cartesian(this_center, theta0, t0);
    theta_to_cartesian(this_center, theta1, t1);
    double phi = (theta1 - theta0)*0.5;
    double r = 1. / cos(phi);
    double thetaq = (theta0+theta1)*0.5;
    theta_to_cartesian(this_center, thetaq, r, q);
    
    chocks.emplace_back( this_center, t0, q );
    chocks.emplace_back( this_center, t1, q );

    // works for here, but seems to cause numeric issues later
    // pull_back_q( poly, si, q );
    // pull_back_t( poly, si, t0 );

    loop.push_back(t0);
    loop.push_back(q);
  }
}

void seg_square_to_loop(const Polygon &poly, size_t si, Loop &loop)
{
  loop.push_back(poly[si].p0);
}

void shave_chocks(const Polygon &poly, Chocks &chocks, Loop &loop)
{
  chocks.clear();
  loop.clear();
  chocks.reserve(poly.size()*4);
  loop.reserve(poly.size()*4);
  for (size_t si = 0; si<poly.size(); ++si)
  {
    switch (poly[si].ctype)
    {
      case CIRCLE:
        seg_circle_to_chock( poly, si, chocks, loop);
        break;
        default:
        seg_square_to_loop( poly, si, loop);
    }
  }
}

void bounding_box(const Loop &loop, Point &box_minp, Point &box_maxp)
{
  if (loop.empty())
  {
    return;
  }
  // get limits of square
  box_minp = loop.front();
  box_maxp = box_minp;
  for (auto &s : loop)
  {
    box_minp.x = std::min( box_minp.x, s.x );
    box_maxp.x = std::max( box_maxp.x, s.x );
    box_minp.y = std::min( box_minp.y, s.y );
    box_maxp.y = std::max( box_maxp.y, s.y );
  }
}

double text_scale(const Point &box_minp, const Point &box_maxp)
{
  double el = norm(box_maxp-box_minp);
  return (1.5 * el) / sqrt(2.);
}

void plot_labels(const Loop &loop, std::ostream &fout)
{
  Point box_minp, box_maxp;
  bounding_box(loop, box_minp, box_maxp);
  auto scale = text_scale(box_minp, box_maxp);

  // loop vertex labels
  for (size_t i=0; i<loop.size(); ++i)
  {
    fout << "\n% loop vertex " << i << "\n";
    fout << "/Times-Bold findfont "<< 0.1*scale << " scalefont setfont\n";
    fout << loop[i].x << " " << loop[i].y << " moveto (" << i << ") show\n";
  }
}
void plot_loop(const Loop &loop, std::ostream &fout, bool do_labels, bool do_fill )
{
  if (loop.empty())
  {
    fout << "\n% loop is empty\n";
    return;
  }
  
  // green shaded filled loop
  if (do_fill)
  {
    fout << "0.25 1 0.25 setrgbcolor\n";
    fout << "\n% loop of " << loop.size() << " points\n";
    fout << "newpath\n";
    fout << loop[0].x << " " << loop[0].y << " moveto\n";
    for (size_t i = 1; i < loop.size(); ++i)
    {
      fout << loop[i].x << " " << loop[i].y << " lineto\n";
    }
    fout << "closepath\n";
    fout << "fill\n\n";
    fout << "0 0 0 setrgbcolor\n";
  }

  // loop
  fout << "\n% loop of " << loop.size() << " points\n";
  fout << "newpath\n";
  fout << loop[0].x << " " << loop[0].y << " moveto\n";
  for (size_t i = 1; i < loop.size(); ++i)
  {
    fout << loop[i].x << " " << loop[i].y << " lineto\n";
  }
  fout << "closepath\n";
  fout << "stroke\n\n";

  if (do_labels)
  {
    plot_labels(loop, fout);
  }
  

}
void print_loop(const Loop &loop)
{
  // loop
  std::cout << "\nloop of " << loop.size() << " points:\n";
  auto old_precision = std::cout.precision(std::numeric_limits<double>::max_digits10);
  for (size_t i = 0; i < loop.size(); ++i)
  {
    std::cout << "  " << loop[i].x << "," << loop[i].y << "\n";
  }
  std::cout.precision(old_precision);
}


inline
void next_uncut(int &i, const std::vector<int> cut_vertices)
{
  if (i >= cut_vertices.size())
  {
    i=0;
  }
  while (cut_vertices[i])
  {
    if (++i >= cut_vertices.size())
    {
      i=0;
    }
  }
}

// true if not a reflex angle at b
bool ccw_triangle( const Point &a, const Point &b, const Point &c)
{
  // vectors
  const Point ba = b-a;
  const Point ca = c-a;
  return clockwise(ca,ba);
}

// true if p is inside triangle abc
bool empty_triangle_strict( const Point &a, const Point &b, const Point &c, const Point &p)
{
  // vectors
  return
  clockwise_strict(b-a,p-a) ||
  clockwise_strict(c-b,p-b) ||
  clockwise_strict(a-c,p-c);
  // if it turns out robustness is an issue for tiny-angle triangles, try also a dot product check.
  // if b-a dot p-a > normsquared(b-a), for two edges of a triangle, then p lies "beyond" the half-spaces trimming the vertices
}

bool empty_triangle( const Point &a, const Point &b, const Point &c, const Point &p)
{
  // vectors
  return
  clockwise(b-a,p-a) ||
  clockwise(c-b,p-b) ||
  clockwise(a-c,p-c);
}

// true if triangle {loop[i], loop[j], loop[k]} does not contain any over loop vertex
bool empty_triangle_strict(int i, int j, int k, const Loop &loop, const std::vector<int> &next_vertex)
{
  auto m = next_vertex[k]; //   efficiency: only need to check non-cut_vertices
  while (m!=i)
  {
    if (!empty_triangle( loop[i], loop[j], loop[k], loop[m] ))
    {
      // if the overlap is m == i or m==k, then consider the triangle empty. But not if m==j the middle vertex
      if (loop[i] != loop[m] || loop[k] != loop[m] )
      {
        return false;
      }
    }
    m = next_vertex[m];
  }
  return true;
}

bool empty_triangle(int i, int j, int k, const Loop &loop, const std::vector<int> &next_vertex)
{
  auto m = next_vertex[k];
  while (m!=i)
  {
    if (!empty_triangle( loop[i], loop[j], loop[k], loop[m] ))
    {
      return false;
    }
    m = next_vertex[m];
  }
  return true;
}


void circumcenter_old(const Point &a, const Point &b, const Point &c, Point &o, double &R )
{
  // https://mathworld.wolfram.com/Circumcenter.html
  // https://mathworld.wolfram.com/Circumradius.html
  
  // formula from
  // https://en.wikipedia.org/wiki/Circumscribed_circle
  const double D = 2. * (a.x * (b.y-c.y) + b.x * (c.y - a.y) + c.x *(a.y-b.y) );
  const double axy2 = (a.x * a.x + a.y * a.y);
  const double bxy2 = (b.x * b.x + b.y * b.y);
  const double cxy2 = (c.x * c.x + c.y * c.y);
  o.x =  (axy2 * (b.y-c.y) + bxy2 * (c.y-a.y) + cxy2 * (a.y-b.y) ) / D;
  o.y =  (axy2 * (c.x-b.x) + bxy2 * (a.x-c.x) + cxy2 * (b.x-a.x) ) / D;
  R = sqrt(((o.x - a.x) * (o.x - a.x) + (o.y - a.y) * (o.y - a.y) +
            (o.x - b.x) * (o.x - b.x) + (o.y - b.y) * (o.y - b.y) +
            (o.x - c.x) * (o.x - c.x) + (o.y - c.y) * (o.y - c.y) ) / 3. );
  // R  =  (abc)/(4rs) for a,b,c side lengths
  // r = (abc)/(4Rs)
  // s = (a+b+c)/2
  // to do: try compute r
}

inline
void circumcenter_segment(const Point &a, const Point &b, Point &o, double &R)
{
  o = (a+b)*0.5;
  const auto na2 = normsquared( o-a );
  const auto nb2 = normsquared( o-b );
  R = sqrt( (na2 + nb2) * 0.5 );
}

void circumcenter(const Point &aa, const Point &bb, const Point &cc, Point &o, double &R )
{
  // check for degenerate triangle
  // new rule
  if (1)
  {
    if (aa==bb || aa==cc)
    {
      circumcenter_segment(bb,cc,o,R);
      return;
    }
    else if (bb==cc)
    {
      circumcenter_segment(aa,bb,o,R);
      return;
    }
  }

  // translate a to origin first
  Point b = bb-aa;
  Point c = cc-aa;
             
  // formula from
  // https://en.wikipedia.org/wiki/Circumscribed_circle
  const double D = 2. * (b.x * c.y + b.y * c.x);
  const double bxy2 = (b.x * b.x + b.y * b.y);
  const double cxy2 = (c.x * c.x + c.y * c.y);
  o.x = aa.x + (  bxy2 * c.y - cxy2 * b.y) / D;
  o.y = aa.y + ( -bxy2 * c.x + cxy2 * b.x) / D;
  R = sqrt(((o.x - b.x) * (o.x - b.x) + (o.y - b.y) * (o.y - b.y) +
            (o.x - c.x) * (o.x - c.x) + (o.y - c.y) * (o.y - c.y) ) / 2. );
}


bool is_empty_circle_b(const Point &a, const Point &b, const Point &c, const Point &m,
                       const Point &o, double R)
{
  if ( normsquared(m-o) > R*R )
  {
    return true;
  }
  if (!clockwise(a-b, m-b) || clockwise(c-b, m-b))
  {
    return true;
  }
  return false;
}


bool is_empty_circle_j(int i, int j, int k, const Loop &loop, const std::vector<int> &next_vertex, const Point &o, double R)
{
  for (int m = 0; m < (int)loop.size(); ++m)
  {
    if (m==i || m==j || m==k)
    {
      continue;
    }
    if (!is_empty_circle_b( loop[i], loop[j], loop[k], loop[m], o, R ))
    {
      return false;
    }
  }
  return true;
}

double triangle_area(const Point &a, const Point &b, const Point &c)
{
  // for accuracy, use the cross produce of the vectors that have the smallest dot product
  // is that a good rule?
  const Point A = c-b;
  const Point B = a-c;
  const Point C = b-a;
  const double AdotB = fabs(dot_product(A,B));
  const double BdotC = fabs(dot_product(B,C));
  const double CdotA = fabs(dot_product(C,A));
  Point cp ; // cross product
  if (AdotB < BdotC)
  {
    if (AdotB < CdotA)
    {
      return 0.5 * fabs(cross_product(A,B));
    }
    // CdotA smallest
    else
    {
      return 0.5 * fabs(cross_product(C,A));
    }
  }
  // not AdotB
  else
  {
    if (BdotC < CdotA)
    {
      return 0.5 * fabs(cross_product(B,C));
    }
    // CdotA smallest
    else
    {
      return 0.5 * fabs(cross_product(C,A));
    }
  }
}

// monotonically increasing function of angle at b
double measure_angle(const Point &a, const Point &b, const Point &c)
{
  const Point A = a-b;
  const Point C = c-b;
  const double Anorm = norm(A);
  const double Cnorm = norm(C);
  const double cross_term = cross_product(C,A) / (Anorm*Cnorm);
  const double dot_term = dot_product(A,C);
  if (dot_term >= 0.)
  {
    if (cross_term >= 0.)
    {
      // angle in [0, pi/2]
      return cross_term;
    }
    else
    {
      // angle in [pi/2, 3pi/2]
      return 4. + cross_term;
    }
  }
  else
  {
    // angle in [3pi/2, 2pi]
    return 2. - cross_term;
  }
}

// circumradius
double circumradius(const Point &pa, const Point &pb, const Point &pc)
{
  // circumradius does not take into account the orientation of the triangle
  const Point A = pc-pb;
  const Point B = pa-pc;
  const Point C = pb-pa;
  const double a = norm(A);
  const double b = norm(B);
  const double c = norm(C);
//  const double s = (a + b + c) / 2.;
//  const double R = (a*b*c) / (4. * sqrt( s * fabs((a+b-s)*(a+c-s)*(b+c-s)) ) );
  const double R = (a*b*c) / sqrt( fabs((a+b+c) * (b+c-a) * (c+a-b) * (a+b-c)) );
  return R;
}

// return the index of a vertex with the smallest angle_measure, and empty of other polygon vertices
//   given k an uncut vertex, traversing via next_vertex
// could switch to a priority queue if this O(k) operation actually takes significant time,
// k = number of vertices, at most 9
int sharpest_empty_triangle(int k, const std::vector<int> &next_vertex,
                            const std::vector<double> &angle_measure,
                            const std::vector<int> &is_empty_triangle,
                            const std::vector<double> &circumradius_measure,
                            const std::vector<int> &empty_circ)
{
  auto j = k;
  // at least two should be empty and non-obtuse
  const double obtuse_thresh = 1.9999;
  const double almost_obtuse_thresh = 1.8;
  bool failed = false;
  while (!is_empty_triangle[j] || angle_measure[j]>=obtuse_thresh || isnan(angle_measure[j]))
  {
    j = next_vertex[j];
    if (j==k) // failsafe
    {
      failed=true;
      break;
    }
  }
  // at least one will always be empty, unless it's totally messed up
  if (failed)
  {
    j=k;
    while (!is_empty_triangle[j])
    {
      j = next_vertex[j];
      if (j==k) // failsafe
      {
        break;
      }
    }
  }
  auto best_j = j;
  int best_empty_triangle = is_empty_triangle[j];
  double best_angle_measure = angle_measure[j];
  double best_circumradius = circumradius_measure[j];
  int best_empty_circ = empty_circ[j];
  bool best_is_bad = !best_empty_triangle || best_angle_measure>=obtuse_thresh || isnan(best_angle_measure);
  j = next_vertex[j];
  const double zero_angle_measure = 1.0e-10;
  do
  {
    // old rule
    //    just desires sharp angle
    //    if (is_empty_triangle[j] && (angle_measure[j] < best_angle_measure))
    
    // old rule, with error: it won't switch off an obtuse angle if the circumradius is good
    // want non-obtuse and small circumradius
    //    if (is_empty_triangle[j] &&
    //        (angle_measure[j] < obtuse_thresh || (angle_measure[j] < best_angle_measure && best_angle_measure >= obtuse_thresh)) &&
    //        (empty_circ[j] > best_empty_circ || (empty_circ[j] == best_empty_circ && circumradius_measure[j] < best_circumradius)))

    // rule
    // 1. empty triangle with non-obtuse angle
    // 1a. zero angle measure, trim off degenerate segments first
    // 2. empty triangle with not-almost-obtuse angle
    // 3a. empty circumcircle
    // 3b. small circumradius
    bool use_j = false;
    // 1.
    const bool j_is_bad = !is_empty_triangle[j] || (angle_measure[j]>=obtuse_thresh) || isnan(angle_measure[j]);
    if (j_is_bad && !best_is_bad)
    {
      // don't use j
    }
    else if (!j_is_bad && best_is_bad)
    {
      // don't use "best"
      // use j regardless of other criteria
      use_j = true;
    }
    else
    {
      // both are bad, or both are good
      assert( j_is_bad == best_is_bad );
      // 0. j has zero angle
      //if (best_angle_measure>0. && angle_measure[j]==0.)
      if (angle_measure[j] <= zero_angle_measure && angle_measure[j] < best_angle_measure)
      {
        use_j = true;
      }
      // 2. best is almost obtuse, and j has a smaller angle
      else if ((isnan(best_angle_measure) && !isnan(angle_measure[j])) ||
           (best_angle_measure>=almost_obtuse_thresh && angle_measure[j]<best_angle_measure ))
      {
        use_j = true;
      }
      // 3. circle criteria
      else if ( empty_circ[j]  > best_empty_circ ||
               (empty_circ[j] == best_empty_circ &&
                (circumradius_measure[j] < best_circumradius ||
                 (isnan(best_circumradius) &&!isnan(circumradius_measure[j])))))
      {
        use_j = true;
      }
    }
    
    if (use_j)
    {
      best_j = j;
      best_empty_triangle = is_empty_triangle[j];
      best_angle_measure = angle_measure[j];
      best_circumradius = circumradius_measure[j];
      best_empty_circ = empty_circ[j];
      best_is_bad = !best_empty_triangle || (best_angle_measure>=obtuse_thresh);
    }
    j = next_vertex[j];
  } while(j!=k);
  return best_j;
}

void update_angle_empty(int j, const Loop &loop,
                        const std::vector<int> &prior_vertex, const std::vector<int> &next_vertex,
                        std::vector<double> &angle_measure, std::vector<int> &is_empty_triangle,
                        std::vector<double> &circumradius_measure, std::vector<int> &empty_circ)
{
  const auto i = prior_vertex[j];
  const auto k = next_vertex[j];
  
  angle_measure[j] = measure_angle( loop[i], loop[j], loop[k] );
  
  // circumradius and center
  // circumradius_measure[j] = circumradius( loop[i], loop[j], loop[k] );
  Point o;
  double R;
  circumcenter( loop[i], loop[j], loop[k], o, R);
  
  empty_circ[j] = is_empty_circle_j(i,j,k, loop, next_vertex, o, R);
  circumradius_measure[j] = R; // to do, replace by R/r
  
  is_empty_triangle[j] = (int) empty_triangle_strict(i,j,k, loop, next_vertex);
}

bool simplify_loop_repeats(Loop &loop)
{
  
  if (loop.size()<2)
  {
    return false;
  }
  bool has_repeat = false;
  {
    Point *o = &loop.back();
    for (auto &p : loop)
    {
      if (p==*o)
      {
        has_repeat=true;
        break;
      }
      o = &p;
    }
  }
  if (has_repeat)
  {
    Loop new_loop;
    new_loop.reserve(loop.size());
    {
      Point *o = &loop.back();
      for (auto &p : loop)
      {
        if (p!=*o)
        {
          new_loop.push_back(p);
          o = &p;
        }
      }
    }
    if (new_loop.empty())
    {
      new_loop.push_back(loop.front());
    }
    new_loop.swap(loop);
    return true;
  }
  return false;
}

// check for polygon having crossings
bool simplify_loop_crossings(Loop &loop)
{
  return false;
}


void simplify_loop(Loop &loop)
{
  // check for polygon having crossings
  simplify_loop_crossings(loop);

  
  // check for repeated vertices
  simplify_loop_repeats(loop);

}

int triangulate_loop( const Loop &loop, Triangles &triangles)
{
  // loops have only a dozen segments, order of magnitude, so efficiency is not an issue here
  triangles.clear();

  int num_uncut = (int) loop.size();
  if (num_uncut<3)
  {
    std::cout << "ERROR: loop of only " << num_uncut << " vertices\n";
    assert(num_uncut>=3);
    return 2;
  }

  triangles.reserve((loop.size()-2)*3);
  
  if (num_uncut==3)
  {
    triangles.push_back(0);
    triangles.push_back(1);
    triangles.push_back(2);
    return 0;
  }
  
  // linked list
  
  // next_vertex = [1...n 0]
  std::vector<int> next_vertex(loop.size());
  std::iota( next_vertex.begin(), next_vertex.end()-1, 1);
  next_vertex.back() = 0;

  // prior_vertex = [n 0..(n-1)]
  std::vector<int> prior_vertex(loop.size());
  std::iota( prior_vertex.begin()+1, prior_vertex.end(), 0);
  prior_vertex.front() = num_uncut-1;

  // angle at each vertex
  // i = index of first vertex of candidate ear triangle, j = 2nd, k = 3rd
  std::vector<double> angle_measure(loop.size());
  std::vector<int> is_empty_triangle(loop.size());
  std::vector<double> circumradius_measure(loop.size());
  std::vector<int> empty_circ(loop.size());
  for (int j = 0; j< (int)loop.size(); ++j)
  {
    update_angle_empty(j, loop, prior_vertex, next_vertex, angle_measure, is_empty_triangle, circumradius_measure, empty_circ);
  }
  
  // verify / debug
  int error_code = 0;
  if (1)
  {
    int num_ok = 0;
    for (int j = 0; j< (int)loop.size(); ++j)
    {
      if (angle_measure[j]<2. && is_empty_triangle[j])
      {
        num_ok++;
      }
    }
    if (num_ok==1)
    {
      std::cout << "Warning: polygon has only one ear to cut." << std::endl;
      error_code = 1;
    }
    if (num_ok==0)
    {
      std::cout << "ERROR: polygon has no ears to cut!!!" << std::endl;
      return 2;
    }
  }
  
  {
    int i, j, k(2);
    
    while (true)
    {
      j = sharpest_empty_triangle(k, next_vertex, angle_measure, is_empty_triangle, circumradius_measure, empty_circ);
      
      // save ear j
      i = prior_vertex[j];
      k = next_vertex[j];
      triangles.push_back(i);
      triangles.push_back(j);
      triangles.push_back(k);
      
      // only one triangle left?
      if (--num_uncut == 3)
      {
        triangles.push_back(k);
        triangles.push_back(next_vertex[k]);
        triangles.push_back(i);

        assert( next_vertex[k] == prior_vertex[i] );
        assert( next_vertex[ next_vertex [k]] == i );
        assert( prior_vertex[ prior_vertex [i]] == k );

        return error_code;
      }
      
      // cut ear j
      next_vertex[i] = k;
      prior_vertex[k] = i;
            
      // update criteria
      update_angle_empty(i, loop, prior_vertex, next_vertex, angle_measure, is_empty_triangle, circumradius_measure, empty_circ);
      update_angle_empty(k, loop, prior_vertex, next_vertex, angle_measure, is_empty_triangle, circumradius_measure, empty_circ);
    }
  }
  // unreachable
  return error_code;
}

int triangulate_loop_simple( const Loop &loop, Triangles &triangles)
{
  triangles.clear();
  
  int num_uncut = (int) loop.size();
  if (num_uncut<3)
  {
    // debugging only, not a real problem
    if (/* DISABLES CODE */ (false))
    {
      std::cout << "ERROR: loop of only " << num_uncut << " vertices\n";
      assert(num_uncut>=3);
      return 2;
    }
    return 0;
  }
  
  triangles.reserve((loop.size()-2)*3);
  
  if (num_uncut==3)
  {
    triangles.push_back(0);
    triangles.push_back(1);
    triangles.push_back(2);
    return 0;
  }
  
  // ersatz linked list
  // next_vertex = [1...n 0]
  std::vector<int> next_vertex(loop.size());
  std::iota( next_vertex.begin(), next_vertex.end()-1, 1);
  next_vertex.back() = 0;
  
  {
    // candidate triangle indices
    int i=0;
    int j=1;
    int k=2;
    
    int max_fails = (int) loop.size();
    while (true)
    {
      const auto cp = cross_product(loop[i]-loop[j], loop[k]-loop[j]);
      if (max_fails==0 || cp==0. || (cp<0. && empty_triangle(i,j,k, loop, next_vertex)))
      {
        // save ear j
        triangles.push_back(i);
        triangles.push_back(j);
        triangles.push_back(k);

        // shortcut if one triangle left
        if (--num_uncut == 3)
        {
          triangles.push_back(k);
          triangles.push_back(next_vertex[k]);
          triangles.push_back(i);
          return 0;
        }

        // cut
        next_vertex[i] = k;
        j = k;
        k = next_vertex[k];
        
        max_fails = (int) loop.size();
      }
      else
      {
        // iterate
        i = j;
        j = k;
        k = next_vertex[k];
        max_fails--;
      }
    }
  }
  return 1; // unreachable
}

void plot_triangles(const Loop &loop, const Triangles &triangles, std::ostream &fout, bool do_labels)
{
  // loop
  fout << "\n% " << (triangles.size()/3) << " triangles of loop of " << loop.size() << " points\n";
  
  // each triangle should have 3 points
  assert(triangles.size() % 3 == 0);
  for (int t = 0; t < triangles.size(); t += 3 )
  {
    fout << "\n% triangle " << (int) (t / 3) << "\n";
    fout << "newpath\n";
    const Point &ta = loop[ triangles[t] ];
    fout << ta.x << " " << ta.y << " moveto\n";
    const Point &tb = loop[ triangles[t+1] ];
    fout << tb.x << " " << tb.y << " lineto\n";
    const Point &tc = loop[ triangles[t+2] ];
    fout << tc.x << " " << tc.y << " lineto\n" <<
    "closepath\n" << "stroke\n\n";
  }

  if (do_labels)
  {
    Point box_minp, box_maxp;
    bounding_box(loop, box_minp, box_maxp);
    auto scale = text_scale(box_minp, box_maxp);
    
    plot_labels(loop, fout);
    // triangle vertex labels
    for (size_t t=0; t<triangles.size(); t+=3)
    {
      const auto ti = t/3;
      fout << "\n% triangle vertex " << ti << "\n";
      fout << "/Times-Bold findfont " << 0.06*scale << " scalefont setfont\n";
      Point barycenter = loop[ triangles[t] ] + loop[ triangles[t+1] ] +  loop[ triangles[t+2] ];
      barycenter.x /= 3;
      barycenter.y /= 3;
      fout << barycenter.x << " " << barycenter.y << " moveto (" << ti << ") show\n";
    }
  }
}


void plot_triangulation(const Loop &loop, const Triangles &triangles, bool do_plot_box, Point bbox_lo, Point bbox_hi, std::string fname,
                        bool do_labels)
{
  fname.append(".ps");
  std::cout << "Writing triangulation to " << fname << std::endl;
  
  std::ofstream fout;
  fout.open (fname);
  
  plot_preamble(fout);
  
  if (do_plot_box)
  {
    plot_bounding_box( bbox_lo, bbox_hi, true, fout);
  }
  
  plot_triangles(loop, triangles, fout);
  
  fout << "\nshowpage\n";
  
  fout.close();
}
                    
