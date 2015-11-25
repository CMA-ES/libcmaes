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
  Matrix3f A;
Vector3f b;
A << 1,2,3,  0,5,6,  0,0,10;
b << 3, 3, 4;
cout << "Here is the matrix A:" << endl << A << endl;
cout << "Here is the vector b:" << endl << b << endl;
Vector3f x = A.triangularView<Upper>().solve(b);
cout << "The solution is:" << endl << x << endl;

  return 0;
}
