//
//  Geometry.cpp
//  AvoidReject
//
//  Created by Mitchell, Scott A on 7/28/21.
//  Copyright © 2021 Mitchell, Scott A. All rights reserved.
//

#include "Geometry.hpp"

#include <assert.h>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>

void chock_to_cartesian( const Point&c, const Point&t, const Point&q, double phi, double r, double &x, double &y )
{
  // add to phi the rotation of t about c
  const Point ct( t-c ); // vector
  // const double scale = norm( ct );  // not needed if Poisson-disk radius is 1
  const double scale = 1.;
  const Point cq( q-c ); // vector
  
  double rot = atan2(ct);
  
  // cannonical orientation is c=(0,0), and t = (1,0), and q is (1,>0)
  // determine cw or ccw orientation of q
  const bool ccw = !clockwise(ct,cq);
  rot += (ccw ? phi : -phi);
  
  x = scale * r * cos(rot) + c.x;
  y = scale * r * sin(rot) + c.y;
}

// intersect two circles, assuming centers are outside the other circle, assume radii = 1
void intersect_circle(const Point &d0, const Point &d1, bool &disjoint, Point &n0, Point &n1)
{
  // assume both radii are 1
  // assert( d0.r == 1. );
  // assert( d1.r == 1. );
  
  Point d01( d1 - d0 ); // vector from d0 center to d1 center
  const double center_distance2 = normsquared(d01); // center distance squared
  assert( center_distance2 > 1. - 1e-4); // we shouldn't be intersecting circles that contain each others' centers
  const double d2 = center_distance2 / 4.; // midpoint_distance squared
  const double h2 = 1. - d2; // height of intersection point above/below midpoint. Uses radius assumption
  if (h2<=0.)
  {
    disjoint = true;
    return;
  }
  disjoint = false;
  
  Point m = (d1 + d0) * 0.5; // midpoint between centers
  Point mPerp( d01.y, -d01.x); // perpedicular to centers vector, of length center_distance
  
  const double perpfactor = sqrt( h2 / (4. * d2) );
  
  n0 = m + mPerp * perpfactor;
  n1 = m - mPerp * perpfactor;

}

// intersect circle with vertical line at x=x
void intersect_x(const Circle &d, double x, bool &disjoint, Point &n0, Point &n1)
{
  const double dx = x - d.c.x;
  const double dy2 = 1. - dx*dx;
  if ( dy2 <= 0. )
  {
    disjoint = true;
    return;
  }
  disjoint = false;
  const double dy = sqrt(dy2);
  n0.x = x;
  n0.y = d.c.y - dy;
  
  n1.x = x;
  n1.y = d.c.y + dy;
}

// intersect circle with horizontal line at y=y
void intersect_y(Circle d, double y, bool &disjoint, Point &n0, Point &n1)
{
  // just swap x-y and call the vertical version, then swap x-y of the resulting points
  std::swap(d.c.x, d.c.y);
  intersect_x( d, y, disjoint, n0, n1 );
  if (!disjoint)
  {
    std::swap( n0.x, n0.y );
    std::swap( n1.x, n1.y );
  }
}

// change a to be in [b,c], by a +/- 2pi
// rotate angle a by 2pi to try to make it lie between b and c, if possible
bool make_between( double &a, double b, double c)
{
  assert(b<=c);
  if (b<=a && a<=c)
  {
    return true;
  }
  double aa=a;
  while (aa>c)
  {
    aa -= 2*M_PI;
  }
  while (aa<b)
  {
    aa += 2*M_PI;
  }
  if (b<=aa && aa<=c)
  {
    a=aa;
    return true;
  }
  return false;
}

// rotate a by 2pi to get it is close to b as possible
void make_close(double &a, double b)
{
  while ( fabs(a - 2*M_PI - b) < fabs(a - b) )
    a -= 2*M_PI;
  while ( fabs(a + 2*M_PI - b) < fabs(a - b) )
    a += 2*M_PI;
}

int scoop_segment(const Segment &seg, const Point &d, Polygon &poly, bool &seg_inside_disk, bool &seg_outside_disk)
{
  seg_inside_disk = false;
  seg_outside_disk = false;
  
  bool disjoint;
  Point n0, n1;
  auto &dp = seg.dp;
  auto &ctype = seg.ctype;
  switch (ctype)
  {
    case NADA:
      return 0;
    case CIRCLE:
      intersect_circle(d, dp, disjoint, n0, n1);
      break;
    case XPLUS:
    case XMINUS:
      intersect_y(d, dp.y, disjoint, n0, n1);
      break;
    case YPLUS:
    case YMINUS:
      intersect_x(d, dp.x, disjoint, n0, n1);
      break;
    default:
      break;
  }
  if (disjoint)
  {
    // return whole seg
    seg_outside_disk = true;
    poly.push_back(seg);
    return 1;
  }
  
  // for circles, convert the points to increasing "angles"
  
  // for lines, order p? by increasing relevant coordinate
  double q0(0.), q1(0.), o0(0.), o1(0.);
  auto &p0 = seg.p0;
  auto &p1 = seg.p1;
  switch (ctype)
  {
    case XPLUS:
    case XMINUS:
      q0 = std::min(p0.x, p1.x);
      q1 = std::max(p0.x, p1.x);
      o0 = n0.x;
      o1 = n1.x;
      break;
    case YPLUS:
    case YMINUS:
      q0 = std::min(p0.y, p1.y);
      q1 = std::max(p0.y, p1.y);
      o0 = n0.y;
      o1 = n1.y;
      break;
    case CIRCLE:
    {
      // segment entirely inside the scooping circle d ?
      {
        const bool p0out = outside(d, p0);
        const bool p1out = outside(d, p1);
        if (!p0out && !p1out)
        {
          // NADA
          seg_inside_disk = true;
          return 0;
        }
      }
      
      // order p0, p1, n0, n1 by angle around segment circle
      double p0angle, p1angle;
      assert( clockwise( dp, p0, p1) );
      angles(dp, p1, p0, p1angle, p0angle);
      
      // p1angle < p0angle because polygon is ccw from interior
      double n0angle = atan2( n0-dp );
      bool n0fixed = make_between(n0angle,p1angle,p0angle);
      double n1angle = atan2( n1-dp );
      bool n1fixed = make_between(n1angle,p1angle,p0angle);
      
      // ensure the difference between the n0 and n1 angles is less than pi
      if (!n1fixed)
      {
        make_close(n1angle,n0angle);
      }
      else if (!n0fixed)
      {
        make_close(n0angle,n1angle);
      }

      // order n0angle < n1angle
      if (n1angle < n0angle)
      {
        std::swap(n0,n1);
        std::swap(n0angle,n1angle);
      }
      
      q0 = p1angle;
      q1 = p0angle;
      o0 = n0angle;
      o1 = n1angle;
    }
      break;
    default:
      break;
  }
  // having exactly one inequality be strict results in not missing intersections
  //    of a disk passing exactly through a segment endpoint
  if (q1<=o0 || q0>=o1) // non-strict
  {
    // disjoint outside
    seg_outside_disk = true;
    poly.push_back(seg);
    return 1;
  }
  // if (q0>=o0 && q1<=o1) // make inequality to purposely generate bad poly cases 
  if (q0>o0 && q1<o1) // strict inequality
  {
    // disjoint inside
    // NADA
    seg_inside_disk = true;
    return 0;
  }

  bool skip_s0(false), skip_s2(false);
  int num_new = 3;
  const bool middle = q0<=o0 && q1>=o1;
  // if middle we retain all three segments
  if (!middle)
  {
    // we only retain 2 segments: s0 s1 or s1 s2
    const bool left = !middle && q0<=o0;
    // right = !middle && !left

    skip_s2 =
    (left  && (ctype == XPLUS  || ctype == YPLUS)) ||
    (!left && (ctype == XMINUS || ctype == YMINUS || ctype == CIRCLE));
    
    skip_s0 = !skip_s2;
    num_new = 2;
  }
  
  // 3 segs: q0n0  n0n1  n1q1
  // where the caller decides if n0n1 is just added or splits the polygon
  // s0.ctype = ctype;
  // s0.dp = dp;
  // s1.ctype = CIRCLE;
  // s1.dp = d;
  // s2.ctype = ctype;
  // s2.dp = dp;
  switch (ctype)
  {
    case XPLUS:
    case YPLUS:
      // Segment(Point q0, Point q1, CurveType dtype, Point d) : p0(q0), p1(q1), ctype(dtype), dp(d) {}

      // s0.p0 = p0;
      // s0.p1 = n0;
      if (!skip_s0)
        poly.emplace_back(p0, n0, ctype, dp);

      // s1.p0 = n0;
      // s1.p1 = n1;
      poly.emplace_back(n0, n1, CIRCLE, d);

      // s2.p0 = n1;
      // s2.p1 = p1;
      if (!skip_s2)
        poly.emplace_back(n1, p1, ctype, dp);
      break;
    case XMINUS:
    case YMINUS:
    case CIRCLE:
      // s0.p0 = p0;
      // s0.p1 = n1;
      if (!skip_s0)
      poly.emplace_back(p0, n1, ctype, dp);
      
      // s1.p0 = n1;
      // s1.p1 = n0;
      poly.emplace_back(n1, n0, CIRCLE, d);

      // s2.p0 = n0;
      // s2.p1 = p1;
      if (!skip_s2)
        poly.emplace_back(n0, p1, ctype, dp);

      break;
    default:
      break;
  }
   
  return num_new;
}

void find_span_end( const size_t span_start, size_t &span_end, const Point &d, const Polygon &pr)
{
  // span start is d, next is not d
  span_end=span_start;
  {
    const Segment *s = &pr[span_end]; // unused value
    do
    {
      span_end++;
      if (span_end>=pr.size())
        span_end = 0;
      s = &pr[span_end];
    } while (! s->is_circle(d) );
  }
}

bool find_span_start(size_t &i, const Point &d, const Polygon &pr)
{
  if (i>=pr.size())
  {
    i = 0;
  }
  
  // advance to first incidence of segment of d
  const Segment *s = &pr[i];
  while (! s->is_circle(d) )
  {
    ++i;
    if (i==pr.size())
    {
      return false;
    }
    else
    {
      s = &pr[i];
    }
  }
  // advance to last consecutive indicence of segment of d
  while ( s->is_circle(d) )
  {
    ++i;
    if (i==pr.size())
    {
      if (pr.front().is_circle(d))
        return false; // we already found the zeroth one earlier
      else
        break; // we haven't found the zeroth one already, starts at the very end pr
    }
    s = &pr[i];
  }
  --i; // i points to last incidence of d before a non-d or the end of pr
  return true;
}

void print_point(const Point &p)
{
  std::cout << "(" << p.x << "," << p.y << ") ";
}

void print_poly(const Polygon &poly)
{
  if (poly.empty())
  {
    std::cout << " polygon is empty" << std::endl;
  }
  else
  {
    std::cout << " polygon has " << poly.size() << " segments." << std::endl;
    auto prec = std::cout.precision();
    std::cout << " precision is " << prec << std::endl;
    for (size_t i=0; i<poly.size(); ++i )
    {
      const auto &seg = poly[i];
      std::cout << "  " << i << ": ";
      switch (seg.ctype)
      {
        case NADA:
          std::cout << "nada ";
          break;
        case CIRCLE:
          std::cout << "circle ";
          std::cout << "center ";
          print_point(seg.dp);
        {
          bool cw = clockwise( seg.dp, seg.p0, seg.p1);
          double angle_start, angle_end;
          // poly is oriented large angle to small angle, clockwise
          if (cw)
            angles( seg.dp, seg.p1, seg.p0, angle_end, angle_start );
          else
            angles( seg.dp, seg.p0, seg.p1, angle_end, angle_start);
          
          // radians to degrees
          angle_start *= 180. * M_1_PI;
          angle_end   *= 180. * M_1_PI;
          std::cout << ", angles (" << angle_start << "," << angle_end << ") ";
        }
          break;
        case XPLUS:
          std::cout << "xplus ";
          break;
        case XMINUS:
          std::cout << "xminus ";
          break;
        case YPLUS:
          std::cout << "yplus ";
          break;
        case YMINUS:
          std::cout << "yminus ";
          break;
      }
      print_point(seg.p0);
      print_point(seg.p1);
      std::cout << std::endl;
    }
  }
}

// check if entire poly is covered by the disk
bool poly_covered(const Polygon &poly, const Point &d)
{
  for (auto &s : poly )
  {
    auto &p = s.p0;
    if (point_outside_circle(p,d))
    {
      return false;
    }
  }
  return true;
}

void scoop_poly(const Polygon &poly, const Point &d, std::vector<Polygon> &result, bool do_debug)
{
  if (poly.empty())
    return;
  
  if (do_debug)
  {
    std::cout << "scoop initial poly:\n";
    print_poly(poly);
  }
  
  result.resize( result.size()+1 );
  result.back().reserve(poly.size()*2);
  // Segment s0, s1, s2;
  
  // keep track of empty intersection types, seg wholely inside or outside the disk
  bool inside(false), outside(false), prior_in(false), prior_out(false), zero_in(false), zero_out(false);
  // for (auto &s : poly)
  for (size_t i=0; i<poly.size(); ++i)
  {
    auto &s = poly[i];

    //s.scoop_segment(d, s0, s1, s2);
    auto num_polys = scoop_segment(s, d, result.back(), inside, outside);
    
    // check for missed topological crossing
    // never encountered in current scoop version 17 Nov 2021, but encountered for test_grid_3(true) periodic
    //   did we go between outside and inside with no crossing?
    if ((inside && prior_out) || (outside && prior_in))
    {
//      s0.ctype = CIRCLE;
//      s0.dp = d.c;
//      s0.p0 = s.p0;
//      s0.p1 = s0.p0;
      result.back().emplace_back( s.p0, s.p0, CIRCLE, d );
      if (outside && prior_in)
      {
        assert( num_polys > 0 );
        assert( result.back().size() > 1 );
        // swap last two segments, as they are out of order
        std::swap( *(result.back().end()-1), *(result.back().end()-2) );
      }
    }
    if (i==0)
    {
      zero_in = inside;
      zero_out = outside;
    }
    prior_in = inside;
    prior_out = outside;
  }
  // check missing intersection between last and first
  // never encountered in current scoop version 17 Nov 2021
  if ((zero_in && outside) || (zero_out && inside))
  {
    auto &s = poly.front();
    result.back().emplace_back( s.p0, s.p0, CIRCLE, d );
    // either case, it is in the right order
  }
  // disk ate the whole polygon
  if (result.back().empty())
  {
    assert(poly_covered(poly,d));
    result.pop_back();
    return;
  }
  // Split result.back() into its connected polygons
  // Either it is the original polygon, with d disjoint, or
  // each connected polygon has a d-segment.
  else
  {
    assert(!poly_covered(poly,d));
    if (do_debug)
    {
      std::cout << "scoop raw result:\n";
      print_poly(result.back());
    }

    size_t new_result_i = result.size()-1;
    
    // Find spans of not-Circle-d-segments
    bool done = false;
    bool first_time = true;
    size_t span_start = 0;
    size_t first_span_start = 0;
    do
    {
      // find start of span
      //   advance to disk, not going past the end of the vector
      if (!find_span_start(span_start, d, result.back()))
      {
        // if there was no span start, then d was not anywhere in result.back()
        // if there was never a span start, then result.back() is unchanged
        if (first_time)
        {
          return;
        }
        // else we just went past the end of the vector
        assert( !result.empty() );
        result.pop_back();
        done = true;
        break;
      }
      else
      {
        if (first_time)
        {
          first_span_start = span_start;
        }
        // find span end
        size_t span_end;
        find_span_end(span_start,span_end,d,result.back());
        // copy start to end to new result vector
        // cases:
        // start == end, this only occurs for a disk scooping a near side of a square, copy the whole thing, done
        // start < end, copy from start to end
        // start > end, copy from start to result.back().end(), and result.back().begin() to end
        assert(span_start>=0);
        assert(span_end>=0);
        assert(span_start<result.back().size());
        assert(span_end<result.back().size());

        if (span_start==span_end)
        {
          if (do_debug)
          {
            std::cout << "scoop span_start==span_end, done\n";
            print_poly(result.back());
          }
          return;
        }
        else if (span_start<span_end)
        {
          if (span_start==0 && span_end==result.back().size()-1)
          {
            // one component
            done = true;
          }
          else if (first_span_start==0 && span_end==result.back().size()-1)
          {
            // last component
            assert(span_start<result.back().size());
            std::rotate( result.back().begin(), result.back().begin()+span_start, result.back().end());
            if (span_start!=span_end+1)
            {
              result.back().resize( span_end - span_start + 1);
            }
            done = true;
          }
          else
          {
            // multiple connected components
            // add subset of result.back() to a different polygon
            result.emplace_back( result.back().begin()+span_start, result.back().begin()+span_end+1 );
            // move old result.back() to the new back
            std::swap( result.back(), *(result.end()-2) );
            // get ready for next search
            span_start = span_end;
          }
        }
        else // span_start > span_end, wrap around end of vector and quit
        {
          std::rotate( result.back().begin(), result.back().begin()+span_start, result.back().end());
          if (span_start!=span_end+1)
          {
            result.back().resize( result.back().size() + (-span_start + span_end + 1) );
          }
          done = true;
        }
      }
      first_time = false;
    }
    while (!done);
    
    // if we got to here, then the result has one or more components
    // each starts with one d seg and ends with one d seg
    //    or does not have d at all
    for (size_t ri = new_result_i; ri < result.size(); ++ri )
    {
      auto &ply = result[ri];
      if (do_debug)
      {
        std::cout << "scoop span result:\n";
        print_poly(ply);
      }

      if (ply.size()<2)
      {
        if (do_debug)
        {
          std::cout << "Warning: poly has only two segments.\n";
        }
        continue; // error, but harmless
      }
      Segment &s = ply.back();
      Segment &t = ply.front();
      if ( t.is_circle(d) && s.is_circle(d) )
      {
        // throw away t.p1 and s.p0, make the segments go from t.p0 to s.p1, for both segments
        // s.p1 = t.p1; // not needed
        t.p0 = s.p0;
        // throw away redundant last segment
        ply.pop_back();
      }
      
      if (do_debug)
      {
        std::cout << "scoop final result:\n";
        print_poly(ply);
      }

    }
  }
}

void scoop_polys(std::vector<Polygon> &polys, const Point &d, std::vector<Polygon> &workspace, bool check_covered)
{
  workspace.clear(); // zzyk
  for (auto &poly : polys)
  {
    if (!check_covered || !poly_covered(poly,d))
    {
      scoop_poly(poly, d, workspace,false);
    }
  }
  workspace.swap(polys);
  // std::swap(workspace,polys);
  // workspace contains old polygons
}

void bounding_box(const Polygon &poly, Point &box_minp, Point &box_maxp)
{
  if (poly.empty())
  {
    return;
  }
  // get limits of square
  box_minp = poly.front().p0;
  box_maxp = box_minp;
  for (auto &s : poly)
  {
    box_minp.x = std::min( box_minp.x, s.p0.x );
    box_maxp.x = std::max( box_maxp.x, s.p0.x );
    box_minp.y = std::min( box_minp.y, s.p0.y );
    box_maxp.y = std::max( box_maxp.y, s.p0.y );
  }
}

void plot_square_biters(Point box_minp, Point box_maxp, const std::vector<Point> &bites, size_t num_bites, std::string fname)
{
  fname.append(".ps");
  std::cout << "Writing square and circles to " << fname << std::endl;
  
  std::ofstream fout;
  fout.open (fname);
  
  fout.precision(std::numeric_limits<double>::max_digits10);
  
  // this preamble has different translate
  // display in ps file
  fout << "%!PS\n"
  << "0.035 0.000 0.000 setrgbcolor\n"
  //  << "72 72 scale\n"
  //  << "4.25 5.5 translate\n"
  << "144 144 scale\n"
  << "1.75 2.75 translate\n"
  << "0.0005 setlinewidth\n";


  // expand box
  double dx = 1.01*(box_maxp.x - box_minp.x);
  double dy = 1.01*(box_maxp.y - box_minp.y);
  
  // expand to 5x5 boxes
  Point minp = box_minp;
  Point maxp = box_maxp;
  minp.x -= 2*dx;
  minp.y -= 2*dy;
  maxp.x += 2*dx;
  maxp.y += 2*dy;
  
  // bounding box is expanded squares
  plot_bounding_box(minp, maxp, false, fout);

  // square
  fout << "\n% square \n";
  fout << "newpath\n" <<
  box_minp.x << " " << box_minp.y << " moveto\n" <<
  box_maxp.x << " " << box_minp.y << " lineto\n" <<
  box_maxp.x << " " << box_maxp.y << " lineto\n" <<
  box_minp.x << " " << box_maxp.y << " lineto\n" <<
  "closepath\n" << "stroke\n\n";
  
  for (size_t b = 0; b<=num_bites; ++b)
  {
    fout << "\n% circle " << b << " of " << num_bites << " \n";
    plot_circle(bites[b], fout);
  }
  
  fout << "\nshowpage\n";
  
  fout.close();
}


void plot_poly( const Polygon &poly, std::string fname )
{
  Point minp, maxp;
  bounding_box(poly, minp, maxp);
  plot_poly( poly, fname, minp, maxp );
}

void plot_preamble( std::ofstream &fout )
{
  fout.precision(std::numeric_limits<double>::max_digits10);
  
  // display in ps file
  fout << "%!PS\n"
  << "0.035 0.000 0.000 setrgbcolor\n"
  //  << "72 72 scale\n"
  //  << "4.25 5.5 translate\n"
  << "144 144 scale\n"
  << "2.125 2.75 translate\n"
  << "0.0005 setlinewidth\n";
}

void plot_bounding_box( Point minp, Point maxp, bool expand_box, std::ofstream &fout)
{
  // expand box
  if (expand_box)
  {
    const double dx = std::max( 0.01, 0.1*(maxp.x - minp.x) );
    const double dy = std::max( 0.01, 0.1*(maxp.y - minp.y) );
    
    minp.x -= dx;
    minp.y -= dy;
    maxp.x += dx;
    maxp.y += dy;
  }
  
  // bounding box
  fout << "\n% bounding box\n";
  fout << "newpath\n" <<
  minp.x << " " << minp.y << " moveto\n" <<
  maxp.x << " " << minp.y << " lineto\n" <<
  maxp.x << " " << maxp.y << " lineto\n" <<
  minp.x << " " << maxp.y << " lineto\n" <<
  "closepath\n" << "stroke\n\n";
}


void plot_poly( const Polygon &poly, std::string fname, Point minp, Point maxp )
{
  fname.append(".ps");
  std::cout << "Writing polygon to " << fname << std::endl;
  
  std::ofstream fout;
  fout.open (fname);
  
  plot_preamble(fout);

  plot_bounding_box( minp, maxp, true, fout);

  if (poly.empty())
  {
    std::cout << " polygon is empty" << std::endl;
  }
  else
  {
    std::cout << " polygon has " << poly.size() << " segments." << std::endl;
    
    // draw each segment of the polygon
    fout << "\n% segments\n";
    for (auto &s : poly)
    {
      s.plot_segment(fout);
    }
  }
  
  fout << "\nshowpage\n";
  
  fout.close();
}

void plot_circle( const Point &c, std::ostream &fout )
{
  double angle_start(0), angle_end(360.);
  
  fout << "newpath\n" <<
  c.x << " " << c.y << " " << 1 << " " << angle_start << " " << angle_end << " arc\n" <<
  "stroke\n";
}

void Circle::plot_circle( std::ostream &fout ) const
{
  ::plot_circle(c,fout);
}

void Segment::plot_segment( std::ostream &fout ) const
{
  if (ctype==CIRCLE)
  {
    bool cw = clockwise( dp, p0, p1);
    double angle_start, angle_end;
    // for ps plotting, want angle_start < angle_end
    if (cw)
      angles( dp, p1, p0, angle_start, angle_end);
    else
      angles( dp, p0, p1, angle_start, angle_end);
    
    // radians to degrees
    angle_start *= 180. * M_1_PI;
    angle_end   *= 180. * M_1_PI;
    
    fout << "newpath\n" <<
    dp.x << " " << dp.y << " " << 1 << " " << angle_start << " " << angle_end << " arc\n" <<
    "stroke\n";
  }
  else
  {
    fout << "newpath\n" <<
    p0.x << " " << p0.y << " moveto\n" <<
    p1.x << " " << p1.y << " lineto\n" <<
    "stroke\n";
  }
  // to do: something different with circle arcs
  
}

void plot_chock(const Point &c, const Point &t, const Point &q, std::ostream &fout, const bool do_center)
{
  // chock center 'c'
  if (do_center)
  {
    fout << "\n% chock center 'c'\n";
    fout << "/Times-Bold findfont 0.1 scalefont setfont\n";
    fout << c.x << " " << c.y << " moveto (c) show\n";
    fout << "newpath\n";
    fout << c.x << " " << c.y << " 0.01 0 360 arc fill\n";
  }
  else
  {
    fout << "\n% skipped chock center 'c'\n";
  }
  
  // draw the chock
  // Point z on the boundary
  Point qc( q-c);
  double qcnorm = norm(qc);
  Point tc( t-c );
  double rchock = norm(tc);
  
  qc *= (rchock / qcnorm);
  
  double angle_t = atan2(tc);
  double angle_q = atan2(qc);
  
  // figure out orientation of arc
  if (angle_q>angle_t && angle_q-angle_t > M_PI)
  angle_q -= 2*M_PI;
  if (angle_q<angle_t && angle_t-angle_q > M_PI)
  angle_q += 2*M_PI;
  
  double angle_start = std::min(angle_q, angle_t) * 180. * M_1_PI;
  double angle_end   = std::max(angle_q, angle_t) * 180. * M_1_PI;
  
  fout << "\n% chock\n";
  fout << "newpath\n" <<
  c.x << " " << c.y << " " << rchock << " " << angle_start << " " << angle_end << " arc\n" <<
  q.x << " " << q.y << " lineto\n" <<
  "closepath\nstroke\n\n";
  
  // q.x << " " << q.y << " lineto\n" <<
  // ax  << " " << ay << " lineto\n" <<
}

void plot_chocks(const Chocks &chocks, std::ostream &fout, const bool do_center)
{
  for (auto &chock : chocks)
  {    
    plot_chock(chock.c,chock.t,chock.q,fout,do_center);
  }
}

void plot_full_circles(const Chocks &chocks, std::ostream &fout)
{
  for (auto &ch: chocks)
  {
    auto &c = ch.c;
    fout << "\n% chock\n";
    fout << "newpath\n" <<
    c.x << " " << c.y << " " << 1 << " " << 0 << " " << 360 << " arc\n" << "stroke\n\n";
  }
}



Polygon square_to_poly(double x0, double y0)
{
  const auto a = 1 / sqrt(2.);
  const auto xa = x0+a;
  const auto ya = y0+a;
  return
  {
    Segment(Point(x0,y0), Point(xa,y0), XPLUS,  Point(x0,y0)),
    Segment(Point(xa,y0), Point(xa,ya), YPLUS,  Point(xa,y0)),
    Segment(Point(xa,ya), Point(x0,ya), XMINUS, Point(xa,ya)),
    Segment(Point(x0,ya), Point(x0,y0), YMINUS, Point(x0,ya))
  };
}
