#include "Statistics.h"
#include <cmath>


Statistics_MeanDev::Statistics_MeanDev()
    : RunningSum(0.0)
{
}
    
void
Statistics_MeanDev::DumpOneResult(double result)
{
    RunningSum += result;
    Collector.push_back(result);
}

void
Statistics_MeanDev::Reset()
{
  Collector.clear();
  RunningSum = 0.0;
}

vector<double> Statistics_MeanDev::GetResultsSoFar() const
{
  vector<double> Results(2);
  Results[0] = RunningSum / Collector.size();

  double SumSquaredDeviations(0.0);
  for (unsigned long i = 0; i < Collector.size(); i++) {
    SumSquaredDeviations += (Results[0] - Collector[i]) * (Results[0] - Collector[i]);
  }

  Results[1] = sqrt(SumSquaredDeviations / (Collector.size() - 1)) / sqrt(Collector.size());

  return Results;
}

void
Statistics::DumpOneResult(double result)
{
  Statistics_MeanDev::DumpOneResult(result);
}

void
Statistics::DumpOneResult(double result1, double result2)
{
  Statistics_MeanDev::DumpOneResult(0.5 * (result1 + result2));
}

