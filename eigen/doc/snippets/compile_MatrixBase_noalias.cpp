#include <Eigen/Eigen>
#include <iostream>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif


using namespace Eigen;
using namespace std;

int main(int, char**)
{
  cout.precision(3);
  Matrix2d a, b, c; a << 1,2,3,4; b << 5,6,7,8;
c.noalias() = a * b; // this computes the product directly to c
cout << c << endl;

  return 0;
}
