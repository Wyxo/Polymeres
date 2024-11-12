#ifndef SIMULATION_H
#define SIMULATION_H

#include "Polymer.h"
#include "Obstacle.h"
#include <ctime>
#include <functional>
class Simulation
{
public:
  
  double T, dt, nu, K, a_interaction, sigma, kappa, kappa_wall, k_string, epsilon, theta, rate_theta, Lx, Ly;
  int N_poly, N_wall, N_obstacle, is_periodic;
  string dir;
  vector<Polymer> crowd;
  vector<Obstacle*> obstacles;
  vector<vector<double>> SelfAngle;
  // Constructors
  Simulation(string dir);

  function<double2(double2, double2)> difference;

  double UpdatePosition(double& Theta, int cpt);
  // Interactions
  void NeighbourInteraction();
  //
  void RunSimulation(int N_Iterations, bool save_optimization);
  void SaveSimulation(const int& step, const double& theta, const double& delta);
  string InitFile();

};



#endif