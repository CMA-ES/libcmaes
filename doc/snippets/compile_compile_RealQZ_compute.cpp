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
  MatrixXf A = MatrixXf::Random(4,4);
MatrixXf B = MatrixXf::Random(4,4);
RealQZ<MatrixXf> qz(4); // preallocate space for 4x4 matrices
qz.compute(A,B);  // A = Q S Z,  B = Q T Z

// print original matrices and result of decomposition
cout << "A:\n" << A << "\n" << "B:\n" << B << "\n";
cout << "S:\n" << qz.matrixS() << "\n" << "T:\n" << qz.matrixT() << "\n";
cout << "Q:\n" << qz.matrixQ() << "\n" << "Z:\n" << qz.matrixZ() << "\n";

// verify precision
cout << "\nErrors:"
  << "\n|A-QSZ|: " << (A-qz.matrixQ()*qz.matrixS()*qz.matrixZ()).norm()
  << ", |B-QTZ|: " << (B-qz.matrixQ()*qz.matrixT()*qz.matrixZ()).norm()
  << "\n|QQ* - I|: " << (qz.matrixQ()*qz.matrixQ().adjoint() - MatrixXf::Identity(4,4)).norm()
  << ", |ZZ* - I|: " << (qz.matrixZ()*qz.matrixZ().adjoint() - MatrixXf::Identity(4,4)).norm()
  << "\n";

  return 0;
}

  return 0;
}
