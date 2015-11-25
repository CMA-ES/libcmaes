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
  Matrix3d m = Matrix3d::Random();
cout << "Here is the matrix m:" << endl << m << endl;
Matrix3d inverse;
bool invertible;
double determinant;
m.computeInverseAndDetWithCheck(inverse,determinant,invertible);
cout << "Its determinant is " << determinant << endl;
if(invertible) {
  cout << "It is invertible, and its inverse is:" << endl << inverse << endl;
}
else {
  cout << "It is not invertible." << endl;
}

  return 0;
}

  return 0;
}

  return 0;
}
