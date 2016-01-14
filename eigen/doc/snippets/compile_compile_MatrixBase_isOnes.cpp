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
  Matrix3d m = Matrix3d::Ones();
m(0,2) += 1e-4;
cout << "Here's the matrix m:" << endl << m << endl;
cout << "m.isOnes() returns: " << m.isOnes() << endl;
cout << "m.isOnes(1e-3) returns: " << m.isOnes(1e-3) << endl;

  return 0;
}

  return 0;
}
