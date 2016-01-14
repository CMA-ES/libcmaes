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
cout << "Here are the coefficients on the 1st super-diagonal and 2nd sub-diagonal of m:" << endl
     << m.diagonal<1>().transpose() << endl
     << m.diagonal<-2>().transpose() << endl;

  return 0;
}
