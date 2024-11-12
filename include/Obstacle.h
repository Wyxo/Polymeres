#ifndef OBSTACLE_H
#define OBSTACLE_H
#include "iostream"
#include "Polymer.h"
#include <functional>
class Obstacle : public Agent
{
public:
  virtual void PolymerInteraction(Polymer&, double const) = 0;
};
#endif  