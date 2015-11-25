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
  MatrixXd m(2,3);
m << 1,2,3,4,5,6;
cout << "here's the 2x3 matrix m:" << endl << m << endl;
cout << "let's resize m to 3x2. This is a conservative resizing because 2*3==3*2." << endl;
m.resize(3,2);
cout << "here's the 3x2 matrix m:" << endl << m << endl;
cout << "now let's resize m to size 2x2. This is NOT a conservative resizing, so it becomes uninitialized:" << endl;
m.resize(2,2);
cout << m << endl;

  return 0;
}

  return 0;
}
