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
  MatrixXi m(3, 3);
m << 1, 2, 3,
     4, 5, 6,
     7, 8, 9;
m = (m.array() >= 5).select(-m, m);
cout << m << endl;

  return 0;
}
