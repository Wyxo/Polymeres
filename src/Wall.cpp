#include "Wall.h"

Wall::Wall(double x_start, double y_start, double x_end, double y_end)
{
  edges.first = double2(x_start, y_start);
  edges.second = double2(x_end, y_end);
  vector = edges.second - edges.first;
  length  = vector%vector; 
}

inline bool CCW(double2& A, double2& B, double2& C)
{
  /*
  Check if three points are listed in a counterclockwise order
  ref: https://bryceboe.com/2006/10/23/line-segment-intersection-algorithm/ 
  */
  return (C.second - A.second) * (B.first - A.first) > (B.second - A.second) * (C.first - A.first);
}

inline bool intersect(double2& A, double2& B, double2& C, double2& D)
{
  /*
  Check if the segment [A,B] intersects the segment [C,D]
  ref: https://bryceboe.com/2006/10/23/line-segment-intersection-algorithm/ 
  */
  return (CCW(A, C, D) != CCW(B, C, D)) && (CCW(A, B, C) != CCW(A, B, D));
}

double2 find_intersection_point(double2 X1, double2 X2, double2 Xg, double2 Xd)
{
  /*
  Find the intersection point of two segments [X1, X2] and [Xg, Xd] knowing that
  they intersect.
  t is the solution to s*X1+(1-s)*X2 = t*Xg + (1-t)*Xd 
  */
  double num = X2.second - Xd.second + (Xd.first-X2.first)*(X1.second-X2.second)/(X1.first-X2.first);
  double denom =Xg.second - Xd.second - (Xg.first - Xd.first)*(X1.second-X2.second)/(X1.first-X2.first);
  double t = num/denom;
  return t*Xg + (1-t)*Xd;
}

void Wall::PolymerInteraction(Polymer& poly, double const kappa_wall)
{
  // artificial new edges to take into account the radius of the pedestrian
  // when we look if the polymers intersect the wall
  double2 extended_first = edges.first - ((0.1+poly.radius)/sqrt(length))*vector;
  double2 extended_second = edges.second + ((0.1+poly.radius)/sqrt(length))*vector;
  for (int idx_atom = 1; idx_atom < poly.N; idx_atom++)
  {
    // Find the closest point on the wall of the atom idx_atom of poly
    double2 dX_FirstEdge = poly.atoms[idx_atom] - edges.first;
    // scalar product between wall vector and the pair difference with the first edge
    double cos = (vector%dX_FirstEdge) / length;
    // if < 0 or > 1 it means that the orthogonal projection is not on the wall and we will use
    // the edge as the closest point
    cos = max(0.0, min(1.0, cos));
    // coordinate of the closest point and difference vector
    double2 nearest = edges.first + cos*vector;
    double2 delta = poly.atoms[idx_atom] - nearest;

    // 0.1 is the radius below which the pedestrian aren't confortable with the wall
    double diff = (poly.radius + 0.1) / (!delta + 1e-16) - 1;
    if (diff > 0)
    {
      poly.forces[idx_atom] += kappa_wall * diff * delta;
    }
    // check if one segment of the polymer intersect the "extended" wall 
    if (intersect(poly.atoms[idx_atom-1], poly.atoms[idx_atom], extended_first, extended_second))
    {
      // intersection point between the wall and the polymer segment
      double2 intersection_point = find_intersection_point(poly.atoms[idx_atom-1], poly.atoms[idx_atom], extended_first, extended_second);
      // distance between the intersection and each edge of the wall to find the closest edge
      double2 dX_firstedge = extended_first - intersection_point;
      double2 dX_secondedge = extended_second - intersection_point;
      double2 dX_closest;
      // Move the atoms toward the closest edge
      if (dX_firstedge%dX_firstedge > dX_secondedge%dX_secondedge)
      {
        dX_closest = dX_secondedge;
      }
      else 
      {
        dX_closest = dX_firstedge;
      }
      poly.forces[idx_atom] += kappa_wall * dX_closest;
      poly.forces[idx_atom-1] += kappa_wall * dX_closest;
    }
  }
}
