#include <iostream>
#include "Simulation.h"
#include "progressbar.hpp"
#include <getopt.h>
using namespace std;


/**
 * @brief Array of long options for command-line argument parsing.
 * 
 * This array defines the long options that can be used with the program.
 * Each entry in the array corresponds to a specific command-line option.
 * 
 * Options:
 * - "save_optimization" (short option: 's'): No argument required. Used to save the optimization results.
 * - "T_max" (short option: 't'): Requires an argument. Specifies the time length of the simulation.
 * - "cpt_max" (short option: 't'): Requires an argument. Specifies the maximum amount of steps during optimization.
 * - "path" (short option: 'p'): Requires an argument. Specifies the file path.
 * 
 * The array is terminated with a sentinel value to indicate the end of options.
 */
const struct option long_options[] = 
{
    {"save_optimization", no_argument, nullptr, 's'},
    {"N_opt", required_argument, nullptr, 'N'},
    {"N_opt_max", required_argument, nullptr, 'n'},
    {"path", required_argument, nullptr, 'p'},
    {"progress_bar", required_argument, nullptr, 'b'},
    {nullptr, 0, nullptr, 0}
};


int main(int argc, char *argv[])
{
  // Default values
  int N_opt = 100;
  int N_opt_max = 10000;
  string dir = ".";
  bool save_optimization = false;
  bool progress_bar = false;

  // Parse command-line arguments
  int option_index = 0;
  int option;
  while ((option = getopt_long(argc, argv, "s:N:n:p:b", long_options, &option_index)) != -1) 
  {
    switch (option) 
    {
      case 's':
        save_optimization = true;
        break;
      case 'N':
        N_opt = atoi(optarg);
        break;
      case 'n':
        N_opt_max = atoi(optarg);
        break;
      case 'p':
        dir = optarg;
        break;
      case 'b':
        progress_bar = true;
        break;
      default:
        break;
    }
  }
  // Run the simulation
  Simulation simu(dir);
  int cpt = 0;
  progressbar bar(N_opt);
  while (cpt < N_opt)
  {
    if (progress_bar){bar.update();}
    simu.RunSimulation(N_opt_max , false);
    for (Polymer& poly : simu.crowd)
    {
      for (int i=0; i < poly.N-1; i++){
        poly.atoms[i] = poly.atoms[i+1];
      }
      poly.atoms[poly.N-1] = poly.atoms[poly.N-2] + 1.4*poly.dt*poly.goal;
    }
    ++cpt;
  }
  return 0;
}

