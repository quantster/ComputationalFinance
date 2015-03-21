////////////////////////////////////////////////////////////////////////////////
//
//  Author: Nitish Kanabar
//  Course: Math 623, Spring 2010
//  Assignment: HW 1
//  Compiler Used: 686-apple-darwin10-g++-4.2.1 (GCC) 4.2.1 (Apple Inc. build 5646) (dot 1)
//
////////////////////////////////////////////////////////////////////////////////

#ifndef _MonteCarloSimulator_H_
#define _MonteCarloSimulator_H_

#include <Vanilla1.h>

class MonteCarloSimulator
{
  public:
    enum RandomWalkType {SingleStep, EulerSpot, EulerLogSpot, Milstein};

    MonteCarloSimulator(double S, double r, double vol, double d,
        unsigned long paths, unsigned long steps)
      : spot(S), riskFreeRate(r), volatility(vol), dividendYield(d),
        numberOfPaths(paths), numberOfSteps(steps) {}

    double ExpectedPayoff(VanillaOption&, RandomWalkType);

    void NumberOfPaths(unsigned long paths) { numberOfPaths = paths; }
    unsigned long NumberOfPaths() { return numberOfPaths; }

    void NumberOfSteps(unsigned long steps) { numberOfSteps = steps; }
    unsigned long NumberOfSteps() { return numberOfSteps; }

  private:
    double spot;
    double volatility;
    double dividendYield;
    double riskFreeRate;

    double deltaT;

    unsigned long numberOfPaths;
    unsigned long numberOfSteps;

    double EulerSpotWalk(const VanillaOption&);
    double EulerLogSpotWalk(const VanillaOption&);
    double MilsteinWalk(const VanillaOption&);
};

#endif
