#include "Pillar.h"

Pillar::Pillar(double x, double y, double radius1)
{
  radius = radius1;
  center.first = x;
  center.second = y;
}

void Pillar::PolymerInteraction(Polymer& poly, double const kappa_wall)
{
  for (int idx_atom = 1; idx_atom < poly.N; idx_atom++)
  {
    double2 dX = distance_periodic(center, poly.atoms[idx_atom], poly.Lx, poly.Ly);
    // 0.1 is the distance within the pillar below which the pedestrian isn't confortable
    double diff = (poly.radius + radius) / (!dX + 1e-16) - 1;
    if (diff > 0)
    {
      poly.forces[idx_atom] += (kappa_wall * diff) * dX;
    }
  }
}