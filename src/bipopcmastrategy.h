
#ifndef BIPOPCMASTRATEGY_H
#define BIPOPCMASTRATEGY_H

#include "ipopcmastrategy.h"
#include <random>

namespace libcmaes
{
  template <class TCovarianceUpdate>
  class BIPOPCMAStrategy : public IPOPCMAStrategy<TCovarianceUpdate>
  {
  public:
    BIPOPCMAStrategy(FitFunc &func,
		     CMAParameters &parameters);
    ~BIPOPCMAStrategy();

    void tell();

    int optimize();

  protected:
    void r1();
    void r2();

  private:
    std::mt19937 _gen;
    std::uniform_real_distribution<> _unif;
    double _lambda_def;
    double _lambda_l;
  };
}

#endif
