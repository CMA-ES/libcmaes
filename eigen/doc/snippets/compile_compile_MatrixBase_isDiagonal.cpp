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
  Matrix3d m = 10000 * Matrix3d::Identity();
m(0,2) = 1;
cout << "Here's the matrix m:" << endl << m << endl;
cout << "m.isDiagonal() returns: " << m.isDiagonal() << endl;
cout << "m.isDiagonal(1e-3) returns: " << m.isDiagonal(1e-3) << endl;


  return 0;
}

  return 0;
}
