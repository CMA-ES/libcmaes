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
  VectorXd v(10);
v.resize(3);
RowVector3d w;
w.resize(3); // this is legal, but has no effect
cout << "v: " << v.rows() << " rows, " << v.cols() << " cols" << endl;
cout << "w: " << w.rows() << " rows, " << w.cols() << " cols" << endl;

  return 0;
}

  return 0;
}
