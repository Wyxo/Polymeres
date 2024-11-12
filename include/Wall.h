#ifndef WALL_H
#define WALL_H
#include "iostream"
#include "Obstacle.h"

class Wall : public Obstacle
{
public:
  void PolymerInteraction(Polymer&, double const);
  std::pair<double2, double2>  edges;
  double2 vector;
  double length;
  Wall(double x_start, double y_start, double x_end, double y_end);
};
#endif

