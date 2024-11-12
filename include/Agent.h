#ifndef AGENT_H
#define AGENT_H
#include <vector>
using namespace std;

// ************* DEFINE NEW TYPES ***************
// 2D pair
typedef pair<double,double> double2;
extern double2 operator+(double2 const a, double2 const b);
extern void operator+=(double2& a, double2 const b);
extern void operator-=(double2& a, double2 const b);
extern double2 operator-(double2 const a, double2 const b);
extern double2 operator*(double2 const a, double2 const b); // double2 * double 2
extern double2 operator*(double const coef, double2 const R); // scalar * double2
extern double2 operator/(double2 const a, double const b); // scalar * double2
extern double operator%(double2 const a, double2 const b); // scalar product
extern double operator!(double2 a); //norm of the vector
extern ostream & operator << (ostream &out, const double2 &c);
double2 distance_periodic(double2 const, double2 const, double const, double const);
// Vector of 2D pair
typedef vector<double2> vector2;
extern vector2 operator+(vector2 const a, vector2 const b);
extern vector2 operator-(vector2 const a, vector2 const b);
extern vector2 operator+(double2 const a, vector2 const b);
extern void operator+=(vector2& a, vector2 const b);
extern void operator-=(vector2& a, vector2 const b);
extern vector2 operator*(double const coef, vector2 const a);
extern double operator!(vector2 const a);
extern double operator%(vector2 const a, vector2 const b);
class Agent
{
public:
};

#endif