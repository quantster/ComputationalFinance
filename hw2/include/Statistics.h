
#ifndef STATISTICS_DEV_H_
#define STATISTICS_DEV_H_

#include <vector>

using namespace std;

class Statistics_MeanDev 
{
  public:
    Statistics_MeanDev();

    virtual void DumpOneResult(double result);
    virtual void Reset();
    virtual std::vector<double> GetResultsSoFar() const;

  private:
    double RunningSum;
    std::vector<double> Collector;
};

class Statistics : public Statistics_MeanDev
{
  public:
    Statistics() {}
    virtual void DumpOneResult(double result);
    virtual void DumpOneResult(double result1, double result2);
};

#endif

