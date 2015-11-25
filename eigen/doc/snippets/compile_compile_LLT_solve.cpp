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
  typedef Matrix<float,Dynamic,2> DataMatrix;
// let's generate some samples on the 3D plane of equation z = 2x+3y (with some noise)
DataMatrix samples = DataMatrix::Random(12,2);
VectorXf elevations = 2*samples.col(0) + 3*samples.col(1) + VectorXf::Random(12)*0.1;
// and let's solve samples * [x y]^T = elevations in least square sense:
Matrix<float,2,1> xy
 = (samples.adjoint() * samples).llt().solve((samples.adjoint()*elevations));
cout << xy << endl;

  return 0;
}

  return 0;
}
