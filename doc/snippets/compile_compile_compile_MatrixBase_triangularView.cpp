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
  Matrix3i m = Matrix3i::Random();
cout << "Here is the matrix m:" << endl << m << endl;
cout << "Here is the upper-triangular matrix extracted from m:" << endl
     << Matrix3i(m.triangularView<Eigen::Upper>()) << endl;
cout << "Here is the strictly-upper-triangular matrix extracted from m:" << endl
     << Matrix3i(m.triangularView<Eigen::StrictlyUpper>()) << endl;
cout << "Here is the unit-lower-triangular matrix extracted from m:" << endl
     << Matrix3i(m.triangularView<Eigen::UnitLower>()) << endl;
// FIXME need to implement output for triangularViews (Bug 885)

  return 0;
}

  return 0;
}

  return 0;
}
