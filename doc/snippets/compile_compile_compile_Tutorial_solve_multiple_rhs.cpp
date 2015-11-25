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
  Matrix3f A(3,3);
A << 1,2,3,  4,5,6,  7,8,10;
Matrix<float,3,2> B;
B << 3,1, 3,1, 4,1;
Matrix<float,3,2> X;
X = A.fullPivLu().solve(B);
cout << "The solution with right-hand side (3,3,4) is:" << endl;
cout << X.col(0) << endl;
cout << "The solution with right-hand side (1,1,1) is:" << endl;
cout << X.col(1) << endl;

  return 0;
}

  return 0;
}

  return 0;
}
