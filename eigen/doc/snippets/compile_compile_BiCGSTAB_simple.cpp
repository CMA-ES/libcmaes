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
    int n = 10000;
  VectorXd x(n), b(n);
  SparseMatrix<double> A(n,n);
  /* ... fill A and b ... */ 
  BiCGSTAB<SparseMatrix<double> > solver;
  solver.compute(A);
  x = solver.solve(b);
  std::cout << "#iterations:     " << solver.iterations() << std::endl;
  std::cout << "estimated error: " << solver.error()      << std::endl;
  /* ... update b ... */
  x = solver.solve(b); // solve again
  return 0;
}

  return 0;
}
