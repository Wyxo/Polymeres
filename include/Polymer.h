#ifndef POLYMER_H
#define POLYMER_H
#include "Agent.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <functional>

class Polymer: public Agent
{
public:
  
  double radius;
  double2 start;
  double2 goal;
  double T;
  double dt;
  double alpha_opt;
  double dt_optimisation;
  double Lx;
  double Ly;
  int Npos;
  int Nneg;
  int N ;
  vector2 atoms;
  vector2 speed;
  vector2 forces;
  Polymer();
  Polymer(double2 s, double2 g, double T, double dt = 0.1, double r = 0.225, double Lx=5, double Ly=5);
  void ElasticResultant(double k_string, const function<double2(double2, double2)>& difference);
  void GroundResultant(double nu, double K, double k_string, const function<double2(double2, double2)>& difference);
  void GroundResultant_periodic(double nu, double K, double k_string);
};

#endif