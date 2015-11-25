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
  typedef Matrix<double,4,Dynamic> Matrix4Xd;
Matrix4Xd M = Matrix4Xd::Random(4,5);
Projective3d P(Matrix4d::Random());
cout << "The matrix M is:" << endl << M << endl << endl;
cout << "M.colwise().hnormalized():" << endl << M.colwise().hnormalized() << endl << endl;
cout << "P*M:" << endl << P*M << endl << endl;
cout << "(P*M).colwise().hnormalized():" << endl << (P*M).colwise().hnormalized() << endl << endl;
  return 0;
}

  return 0;
}

  return 0;
}
