
#include <string>
#include <unordered_map>
#include "Polymer.h"
#include "Wall.h"
#include "Pillar.h"
std::unordered_map<std::string, double> parseParameters(const std::string& filename);
std::vector<Polymer> readAgentsFromFile(const std::string& filename, double T, double dt, double Lx, double Ly);
void readEnvironmentFromFile(const std::string& filename, std::vector <Obstacle*>& Obstacles);