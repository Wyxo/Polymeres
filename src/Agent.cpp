#include <iostream>
#include "Agent.h"
#include <cmath>
using namespace std;

// ***********  OPERATORS ON 2D PAIRS  ***********
double2 operator+(double2 const a, double2 const b) {    return double2(a.first+b.first, a.second+b.second);   }
void operator+=(double2& a, double2 const b) 
{
    a.first += b.first; 
    a.second += b.second;
}
double2 operator-(double2 const a, double2 const b) {    return double2(a.first-b.first, a.second-b.second);    }
void operator-=(double2& a, double2 const b) 
{
    a.first -= b.first;
    a.second -= b.second;
}
double2 operator*(double2 const a, double2 const b) {    return double2(a.first*b.first, a.second*b.second);    } // double2 * double 2
double2 operator*(double const coef, double2 const R) {    return double2(coef * R.first, coef * R.second);    } // scalar * double2
double2 operator/(double2 const a, double const b) {    return double2(a.first/b, a.second/b);    }
double operator%(double2 const a, double2 const b) {    return double(a.first*b.first + a.second*b.second);    } // scalar product
double operator!(double2 const a) {  return sqrt(a%a);    } // norm of the vector
ostream & operator << (ostream &out, const double2 &c) {   out << "(" << c.first << "," << c.second << ")";   return out;}
double2 distance_periodic(double2 const a, double2 const b, double const Lx, double const Ly)
{
  double2 dX = b - a;
  dX.first -= round(dX.first/Lx)*Lx;
  dX.second -= round(dX.second/Ly)*Ly;
  return dX;
}
// ***********  OPERATORS ON VECTORS OF 2D PAIRS  ***********
vector2 operator+(vector2 const a, vector2 const b)
{
  vector2 result;
  
  for (int i=0; i<a.size(); i++){result.push_back(a[i]+b[i]);}
  return result;
}
vector2 operator-(vector2 const a, vector2 const b)
{
  vector2 result;
  for (int i=0; i<a.size(); i++){result.push_back(a[i]-b[i]);}
  return result;
}
void operator+=(vector2& a, vector2 const b){for (int i=0; i<a.size(); i++){a[i] += b[i];}
}
void operator-=(vector2& a, vector2 const b){for (int i=0; i<a.size(); i++){a[i] -= b[i];}}

vector2 operator+(double2 const a, vector2 const b)
{
  vector2 result;
  for (int i=0; i<b.size(); i++){result.push_back(b[i]+a);}
  return result;
}
vector2 operator*(double const coef, vector2 const a)
{
  vector2 result;
  result.reserve(a.size());  // Préallouer la mémoire
  for (auto& p : a) {
      result.push_back({coef * p.first, coef * p.second});
  }
  return std::move(result);
}

/*vector2 result;
  for (int i=0; i<a.size(); i++){result.push_back(coef*a[i]);}
  return result;
}*/

double operator!(vector2 const a)
{
  double result = 0;
  for (int i=0; i<a.size(); i++){result += !a[i];}
  return result;
}

double operator%(vector2 const a, vector2 const b)
{
  double result = 0;
  for (int i=0; i<a.size(); i++){result += a[i]%b[i];}
  return result;
}