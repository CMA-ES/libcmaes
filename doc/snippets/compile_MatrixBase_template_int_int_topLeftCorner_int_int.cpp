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
  Matrix4i m = Matrix4i::Random();
cout << "Here is the matrix m:" << endl << m << endl;
cout << "Here is m.topLeftCorner<2,Dynamic>(2,2):" << endl;
cout << m.topLeftCorner<2,Dynamic>(2,2) << endl;
m.topLeftCorner<2,Dynamic>(2,2).setZero();
cout << "Now the matrix m is:" << endl << m << endl;

  return 0;
}
