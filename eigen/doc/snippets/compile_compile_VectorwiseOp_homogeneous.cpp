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
  typedef Matrix<double,3,Dynamic> Matrix3Xd;
Matrix3Xd M = Matrix3Xd::Random(3,5);
Projective3d P(Matrix4d::Random());
cout << "The matrix M is:" << endl << M << endl << endl;
cout << "M.colwise().homogeneous():" << endl << M.colwise().homogeneous() << endl << endl;
cout << "P * M.colwise().homogeneous():" << endl << P * M.colwise().homogeneous() << endl << endl;
cout << "P * M.colwise().homogeneous().hnormalized(): " << endl << (P * M.colwise().homogeneous()).colwise().hnormalized() << endl << endl;
  return 0;
}

  return 0;
}
