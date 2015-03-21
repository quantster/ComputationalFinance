/*
  MonteCarloSimulator
*/

#ifndef MONTECARLO_SIMULATOR_H_
#define MONTECARLO_SIMULATOR_H_

#include "Random3.h"
#include "Option.h"
#include "Statistics.h"

class Simulator
{
  public:
    enum VariatesType {Normal, Antithetic, Control};

    Simulator(unsigned long NumberOfPaths_);

    virtual ~Simulator();

    virtual void RunSimulation(
          UpAndOutCall& option,
          double Spot,
          double Vol,
          double r,
          double d,
          unsigned long NumberOfDates,
          Statistics& gatherer,
          VariatesType);


  protected:
    void SimulateWithNormalVariates(UpAndOutCall& option, 
                  double Spot, double pre0, double pre1, double discount, 
                  MJArray& variate, Statistics& gatherer);

    void SimulateWithAntitheticVariates(UpAndOutCall& option, 
                  double Spot, double pre0, double pre1, double discount, 
                  MJArray& variate, Statistics& gatherer);

    void SimulateWithControlVariates(UpAndOutCall& option,
                  double Spot, double pre0, double pre1, double discount,
                  MJArray& variate, Statistics& gatherer);

  private:
    unsigned long NumberOfPaths;
    RandomRand* Generator;

};

#endif
