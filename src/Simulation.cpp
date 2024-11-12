#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <iomanip>

#include <sys/stat.h>
#include <sys/types.h>
#include "Simulation.h"
#include "file_parser.h"

using namespace std;

double2 euclidian_distance(double2 a, double2 b)
{
  return b - a;
}

double2 periodic_distance(double2 a, double2 b, double Lx, double Ly)
{
  double2 d = b - a;
  d.first -= round(d.first / Lx) * Lx;
  d.second -= round(d.second / Ly) * Ly;
  return d;
}

Simulation::Simulation(string dir1)
{
  dir = dir1;
  unordered_map<string, double> dic = parseParameters(dir);
  T = dic["T"];  // Time window
  dt = dic["dt"];// Time step between two positions
  nu = dic["nu"];// width of the floor field around the goal
  K = dic["K"];  // intensity of the desir to reach the target
  a_interaction = dic["a_interaction"]; // intensity of pedestrian midrange repulsion
  sigma = dic["sigma"]; // length scale of the pedestrian repulsion
  kappa = dic["kappa"]; // intensity of contact repulsion
  kappa_wall = dic["kappa_wall"]; // intensity of repulsion with walls
  k_string = dic["k_string"]; // velocity control
  epsilon = dic["epsilon"]; // convergence threshold
  theta = dic["theta"];     // starting temperature
  rate_theta = dic["rate_theta"]; // rate of temperature  decay
  is_periodic = dic["periodic"];
  Lx = dic["Lx"]; // size of system in case of periodic boundaries
  Ly = dic["Ly"];

  if (is_periodic)
  {
    difference = [](double2 a, double2 b){return euclidian_distance(a, b);};
  }
  else
  {
    difference = [this](double2 a, double2 b) { return periodic_distance(a, b, this->Lx, this->Ly); };
  }
  N_poly = 0;
  // Get time to save data in corresponding directory
  // Load the vector of pedestrian
  crowd = readAgentsFromFile(dir, T, dt, Lx, Ly);
  N_poly = crowd.size();

  // Load the obstacles
  readEnvironmentFromFile(dir, obstacles);
  N_obstacle = obstacles.size();
  SelfAngle.resize(N_poly);
  for (int idx_poly = 0; idx_poly < N_poly; idx_poly++)
  {
    SelfAngle[idx_poly].resize(crowd[idx_poly].N);
  }
}

random_device rd;
mt19937 gen(rd());
normal_distribution<double> d(0, 1);

double Simulation::UpdatePosition(double& Theta, int cpt)
{
  this->Simulation::NeighbourInteraction();
  // max of the gradient of the cost to control the convergence
  double max_F = 0;
  for (Polymer& poly : crowd)
  {
    // FORCES CALCULATION
    /*
    if (is_periodic){cout << 1 << endl;poly.GroundResultant_periodic(nu, K, k_string);}
    else{poly.GroundResultant(nu, K, k_string, difference);}
    */
    poly.ElasticResultant(k_string, difference);
    for (int idx_obstacle = 0; idx_obstacle < N_obstacle; ++idx_obstacle)
    {
      obstacles[idx_obstacle]->PolymerInteraction(poly, kappa_wall);
    }
    // FIRE ALGORITHM
    double P = 0;
    P = poly.forces%poly.speed;
    if (P>0)
    {
      ++poly.Npos;
      poly.Nneg = 0;
      if (poly.Npos > 20)
      {
        poly.dt_optimisation = min(1.1*poly.dt_optimisation, 1e-1);
        poly.alpha_opt *= 0.99;
      }
    }
    else
    {
      ++poly.Nneg;
      poly.Npos = 0;
      if (cpt > 20)
      {
        poly.dt_optimisation = max(1e-7, 0.5*poly.dt_optimisation);
        poly.alpha_opt = 0.25; 
      }
      poly.atoms -= 0.5*poly.dt_optimisation*poly.speed;
      poly.speed = vector2(poly.N, double2(0,0));
    } 
    // SEMI IMPLICIT EULER INTEGRATION
    double normF = !poly.forces;
    max_F += normF/poly.N;
    poly.speed += poly.dt_optimisation*poly.forces;
    if (normF > 0){
      poly.speed = (1-poly.alpha_opt)*poly.speed + poly.alpha_opt*(!poly.speed/normF)*poly.forces;
    }

    poly.atoms += poly.dt_optimisation*poly.speed;
    for (int idx_atom=1; idx_atom < poly.N; idx_atom++)
    {
      poly.atoms[idx_atom] += poly.dt_optimisation*sqrt(Theta)*double2(d(gen), d(gen));
    }
    poly.forces = vector2(poly.N, double2(0,0));
  }
  return max_F/(crowd.size());
}

void Simulation::RunSimulation(int N_Iterations, bool save_optimization)
{
  double delta = 1000;
  int cpt = 0;
  double theta_temp = theta;
  while (delta > epsilon && cpt < N_Iterations)
  {
    // delta represents the maximum of the cost's gradient and
    // caracterize the progress of the convergence
    delta = this -> UpdatePosition(theta_temp, cpt);
    cpt++;

    // Decrease temperature at rate_theta
    if (theta_temp > rate_theta){theta_temp -= rate_theta;    }
    else {theta_temp = 0;}
    if (save_optimization && cpt % 50 == 0){this -> SaveSimulation(cpt, theta_temp, delta);    }
  }
  SaveSimulation(cpt, theta_temp, delta);
  //cout << cpt << endl;
}

// Vision function
inline double W(double theta){return 1/(1 + theta * theta);}

inline double modulo_pi(double angle)
{ 
  if (angle > M_PI)  {return angle - 2*M_PI;}
  else if(angle < -M_PI){return angle + 2*M_PI;}
  return angle;
}

void Simulation::NeighbourInteraction()
{
  // Calculate the local direction of each pedestrian
  for (int idx_poly = 0; idx_poly < N_poly; idx_poly++)
  {
    for (int idx_atom = 1; idx_atom < crowd[idx_poly].N-1; idx_atom++)
    {
      double2 dX = difference(crowd[idx_poly].atoms[idx_atom-1], crowd[idx_poly].atoms[idx_atom]);
      SelfAngle[idx_poly][idx_atom] = atan2(dX.second, dX.first);
    }
  }
  // for each pedestrian
  for (int idx_poly1 = 0; idx_poly1 < N_poly - 1; idx_poly1++)
  {
    // for each neighbour
    for (int idx_poly2 = idx_poly1 + 1; idx_poly2 < N_poly; idx_poly2++)
    {
      // for each atom
      for (int idx_atom = 1; idx_atom < crowd[idx_poly1].N; idx_atom++)
      {
        // pair difference vector
        double2 dX = difference(crowd[idx_poly2].atoms[idx_atom], crowd[idx_poly1].atoms[idx_atom]);
        double dXnorm = !dX + 1e-16;

        // norm of the interaction force without vision taken into account
        double2 F = (a_interaction * exp(-(dXnorm - crowd[idx_poly1].radius - crowd[idx_poly2].radius) / sigma) / dXnorm) * dX;

        // delta angle compared to horizontal 
        /*
        double dAngle1 = atan2(-dX.second, -dX.first);
        double dAngle2 = atan2(dX.second, dX.first);

        // angle at which pedestrian 2 is seen by pedestrian 1 (between -pi and pi)
        double Angle1 = modulo_pi(SelfAngle[idx_poly1][idx_atom] - dAngle1);
        // angle at which pedestrian 1 is seen by pedestrian 2 
        double Angle2 = modulo_pi(SelfAngle[idx_poly2][idx_atom] - dAngle2);

        // weight the force with the visual prefactor
        double W1 = W(Angle1);
        double W2 = W(Angle2);
        */
        crowd[idx_poly1].forces[idx_atom] += F;
        crowd[idx_poly2].forces[idx_atom] -= F;
      
        // contact force taken into account if ri + rj > dPositionNorm
        double diff = (crowd[idx_poly1].radius + crowd[idx_poly2].radius) / dXnorm - 1;
        if (diff>0)
        {
          double2 ContactForce = kappa * diff * dX;
          crowd[idx_poly1].forces[idx_atom] += ContactForce;
          crowd[idx_poly2].forces[idx_atom] -= ContactForce;
        }
      }      
    }
  }
}

string Simulation::InitFile()
{
  auto t = time(nullptr);
  auto tm = *localtime(&t);
  ostringstream oss;
  ostringstream oss2;
  ostringstream oss3;
  oss << put_time(&tm, "%Y-%m-%d %H-%M");
  oss2 << put_time(&tm, "%Y-%m");
  oss3 << put_time(&tm, "%Y-%m-%d");
  string full_date_str = oss.str();
  string month_str = oss2.str();
  string date_str = oss3.str();
  month_str = "saves/" + month_str;
  date_str = month_str + "/" + date_str;
  full_date_str = date_str + "/" + full_date_str;
  mkdir(month_str.c_str(), 0777);
  mkdir(date_str.c_str(), 0777);
  mkdir(full_date_str.c_str(), 0777);
  // Copy each file in the directory
  vector<string> files = {"parameters.txt", "environment.txt", "crowd.txt"};
  for (string filename:files)
  {
    ifstream src(filename.c_str(), ios::binary);
    string dest1 = full_date_str + "/" + filename;
    ofstream dest(dest1.c_str(), ios::binary);
    dest << src.rdbuf();
  }
  string filename = full_date_str + "/" + "data.txt";
  ofstream outFile(filename);
  outFile.close();
return filename;
}

void Simulation::SaveSimulation(const int& step, const double& theta, const double& delta)
{ 
  ofstream out;
  string filename = dir + "/data.txt";
  out.open(filename, ios::app);
  if (!out)
  {
      cerr << "Error while opening: " << filename << endl;
      return;
  }
  out << "ITEM : STEP" << endl;
  out << step << endl;
  out << "ITEM : DELTA" << endl;
  out << delta << endl;
  out << "ITEM : TEMPERATURE" << endl;
  out << theta << endl;
  out << "ITEM : NUMBER OF PEDESTRIAN" << endl;
  out << N_poly << endl;
  out << "ITEM : TIME LENGTH" << endl;
  out << crowd[0].N << endl;
  for (int i = 0;  i < N_poly; i++)
  {
    out << "ITEM : PEDESTRIAN " << i << " " << crowd[i].radius << endl;
    for (int j = 0; j < crowd[i].N; j++)
    {
      out << j*dt << " " <<  crowd[i].atoms[j].first << " " << crowd[i].atoms[j].second << endl;
    }
    out << endl;  
  }
  out << endl;  
  out << endl;  
  out.close();
}
