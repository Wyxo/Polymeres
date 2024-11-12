#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include "file_parser.h"
#include "Polymer.h"
#include "Pillar.h"
#include "Wall.h"

using namespace std;


void readEnvironmentFromFile(const string& dir, vector <Obstacle*>& obstacles)
{
  string filename = dir + "/environment.txt";
  ifstream file(filename);
  if (!file.is_open())
  {
    cerr << "Error while opening the file" << filename << endl;
    return ;
  }
  string line;
  while (getline(file, line))
  {
    if (line.find("PILLAR") != string::npos)
    {
      double x, y, r;
      double temp;
      stringstream ss(line);
      ss >> temp;

      // get the radius
      file >> r;
      // get the coordinates
      file >> x >> y;
      obstacles.push_back(new Pillar(x, y, r));
    }
    if (line.find("WALL") != string::npos)
    {
      stringstream ss(line);
      double temp;
      ss >> temp;
      double x1, y1, x2, y2;
      // get the coordinates
      file >> x1 >> y1;
      file >> x2 >> y2;
      obstacles.push_back(new Wall(x1, y1, x2, y2));
    }
  }
}

vector<Polymer> readAgentsFromFile(const string& dir, double T, double dt, double Lx, double Ly)
{
  string filename = dir + "/crowd.txt";
  vector<Polymer> Polymers;
  ifstream file(filename);

  if (!file.is_open())
  {
    cerr << "Error while opening the file" << filename << endl;
    return Polymers;
  }

  string line;
  while (getline(file, line))
  {
    if (line.find("AGENT") != string::npos)
    {
      Polymer poly;
      stringstream ss(line);
      string temp;
      double2 s;
      double2 g;
      double r;
      // agent ID
      ss >> temp;
      // agent radius
      file >> r;

      // start and goal coordinates
      file >> s.first >> s.second;
      file >> g.first >> g.second;

      Polymers.push_back(Polymer(s, g, T, dt, r, Lx, Ly));
    }
  }
  file.close();
  return Polymers;
}

unordered_map<string, double> parseParameters(const string& dir)
{
  string filename = dir + "/parameters.txt";
  unordered_map<string, double> parameters;
  ifstream file(filename);

  if (!file.is_open())
  {
    cerr << "Error while opening the file" << filename << endl;
    return parameters;
  }

  string line;
  while (getline(file, line))
  {
    // Ignorer les lignes vides ou commentées
    if (line.empty() || line[0] == '#')
    {
      continue;
    }

    istringstream iss(line);
    string key;
    double value;

    // Exctract the key and the value
    if (!(iss >> key >> value))
    {
      cerr << "Conversion error : cannot read the line correcly" << endl;
      continue;
    }

    // Store the parameter
    parameters[key] = value;
  }

  file.close();
  return parameters;
}
