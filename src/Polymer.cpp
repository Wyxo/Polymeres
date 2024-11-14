#include <iostream>
#include <cmath>
#include "Polymer.h"
#include "Agent.h"
using namespace std;

Polymer::Polymer()
: radius(0.225), start{0, 0}, goal{2, 0}, T(5), dt(0.1), N(50)
{
}

Polymer::Polymer(double2 s, double2 g, double T1, double dt1, double r, double Lx1, double Ly1)
{
  radius = r;
  T = T1;
  dt = dt1;
  N = round(T / dt);
  Lx = Lx1;
  Ly = Ly1;
  start = s;
  goal = g;
  forces.reserve(N);
  // optimisation parameters
  alpha_opt = 0.25;
  dt_optimisation = 1e-4;
  Npos = 0;
  Nneg = 0;
  // 1.4 is considered the typical velocity for pedestrian
  speed = vector2(N, double2(0,0));
  forces = vector2(N, double2(0,0));
  atoms.push_back(s);
  double Tmax = !(g-s) / 1.4;
  for (int idx_atom = 1; idx_atom < N; idx_atom++)
  {
    if (idx_atom * dt < Tmax){atoms.push_back(atoms[idx_atom-1] + (dt/Tmax)*(g - s));}
    // Stay still if the goal is reached
    else{atoms.push_back(g);}
  }
}

void Polymer::ElasticResultant(double k_string, const function<double2(double2, double2)>& difference)
{
  for (int idx_atom = 1; idx_atom < N-1; idx_atom++)
  {
    forces[idx_atom] += (k_string / dt) * (difference(atoms[idx_atom], atoms[idx_atom-1]) + difference(atoms[idx_atom], atoms[idx_atom+1]));
  }
  // the last atom is attached only from one side
  forces[N-1] += (k_string / dt) * difference(atoms[N-1],atoms[N-2]);
}

void Polymer::GroundResultant(double nu, double K, double k_string,const function<double2(double2, double2)>& difference)
{  
  double2 dX_goal;
  double dX_norm;
  for (int idx_atom = 1; idx_atom < N-1; idx_atom++)
  {
    dX_goal = difference(atoms[idx_atom], goal);
    dX_norm = !dX_goal;
    forces[idx_atom] += K * nu * exp(-nu * dX_norm)*dX_goal;
  }
  dX_goal = difference(atoms[N-1], goal);
  forces[N-1] += (k_string*1.4+K/1.4)*dX_goal/sqrt(dX_goal%dX_goal+1);
}

void Polymer::GroundResultant_periodic(double nu, double K, double k_string)
{  
  forces[N-1] += (k_string*1.4+K/1.4)*goal;
}