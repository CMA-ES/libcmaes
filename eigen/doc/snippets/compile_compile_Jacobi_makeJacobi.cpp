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
  Matrix2f m = Matrix2f::Random();
m = (m + m.adjoint()).eval();
JacobiRotation<float> J;
J.makeJacobi(m, 0, 1);
cout << "Here is the matrix m:" << endl << m << endl;
m.applyOnTheLeft(0, 1, J.adjoint());
m.applyOnTheRight(0, 1, J);
cout << "Here is the matrix J' * m * J:" << endl << m << endl;
  return 0;
}

  return 0;
}
