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
  std::string sep = "\n----------------------------------------\n";
Matrix3d m1;
m1 << 1.111111, 2, 3.33333, 4, 5, 6, 7, 8.888888, 9;

IOFormat CommaInitFmt(StreamPrecision, DontAlignCols, ", ", ", ", "", "", " << ", ";");
IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
IOFormat OctaveFmt(StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");
IOFormat HeavyFmt(FullPrecision, 0, ", ", ";\n", "[", "]", "[", "]");

std::cout << m1 << sep;
std::cout << m1.format(CommaInitFmt) << sep;
std::cout << m1.format(CleanFmt) << sep;
std::cout << m1.format(OctaveFmt) << sep;
std::cout << m1.format(HeavyFmt) << sep;

  return 0;
}
