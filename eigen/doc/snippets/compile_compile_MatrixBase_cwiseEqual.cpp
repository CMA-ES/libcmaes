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
  MatrixXi m(2,2);
m << 1, 0,
     1, 1;
cout << "Comparing m with identity matrix:" << endl;
cout << m.cwiseEqual(MatrixXi::Identity(2,2)) << endl;
int count = m.cwiseEqual(MatrixXi::Identity(2,2)).count();
cout << "Number of coefficients that are equal: " << count << endl;

  return 0;
}

  return 0;
}
