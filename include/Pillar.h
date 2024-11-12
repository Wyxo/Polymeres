#ifndef PILLAR_H
#define PILLAR_H
#include "iostream"
#include "Obstacle.h"

class Pillar : public Obstacle
{
public:
  void PolymerInteraction(Polymer&, double const);
  double2 center;
  double radius;
  Pillar(double x, double y, double radius);
};
#endif