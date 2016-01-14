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
  int array[24];
for(int i = 0; i < 24; ++i) array[i] = i;
cout << Map<MatrixXi, 0, Stride<Dynamic,2> >
         (array, 3, 3, Stride<Dynamic,2>(8, 2))
     << endl;

  return 0;
}

  return 0;
}
