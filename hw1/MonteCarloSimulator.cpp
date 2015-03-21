////////////////////////////////////////////////////////////////////////////////
//
//  Author: Nitish Kanabar
//  Course: Math 623, Spring 2010
//  Assignment: HW 1
//  Compiler Used: 686-apple-darwin10-g++-4.2.1 (GCC) 4.2.1 (Apple Inc. build 5646) (dot 1)
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include "MonteCarloSimulator.h"
#include "Random1.h"

using namespace std;

//------------------------------------------------------------------------------
// This function performs the simulation and returns the expected payoff of
// the specified option
//------------------------------------------------------------------------------
double 
MonteCarloSimulator::ExpectedPayoff(VanillaOption& option, RandomWalkType walk)
{
  // precalculation
  deltaT = option.GetExpiry() / numberOfSteps;

  double price;
  switch (walk) {
    case SingleStep:
      {
        double steps = numberOfSteps;
        numberOfSteps = 1;
        deltaT = option.GetExpiry();
        price = EulerLogSpotWalk(option);
        numberOfSteps = steps;
        deltaT = option.GetExpiry() / numberOfSteps;
      }
      break;
    case EulerSpot:
      price = EulerSpotWalk(option);
      break;
    case EulerLogSpot:
      price = EulerLogSpotWalk(option);
      break;
    case Milstein:
      price = MilsteinWalk(option);
      break;
    default:
      break;
  };

  return exp(-riskFreeRate * option.GetExpiry()) * price;
}

//------------------------------------------------------------------------------
// This function returns the payoff computed using Euler scheme for spot
//------------------------------------------------------------------------------
double
MonteCarloSimulator::EulerSpotWalk(const VanillaOption& option)
{
  double sumPayoff = 0;

  double pre0 = 1 + (riskFreeRate - dividendYield) * deltaT;
  double pre1 = volatility * sqrt(deltaT);

  for (int i = 0; i < numberOfPaths; i++) {
    double movedSpot = spot;
    for (unsigned long j = 0; j < numberOfSteps; j++) {
      movedSpot = movedSpot * (pre0 + pre1 * GetOneGaussianByBoxMuller());
    }
    sumPayoff += option.OptionPayOff(movedSpot);
  }
  return (sumPayoff / numberOfPaths);
}

//------------------------------------------------------------------------------
// This function returns the payoff using Euler scheme for log spot
//------------------------------------------------------------------------------
double
MonteCarloSimulator::EulerLogSpotWalk(const VanillaOption& option)
{
  double sumPayoff = 0;

  double pre0 = (riskFreeRate - dividendYield - 0.5 * volatility * volatility) * deltaT;
  double pre1 = volatility * sqrt(deltaT);

  for (unsigned long i = 0; i < numberOfPaths; i++) {
    double movedSpot = spot; 
    for (unsigned long j = 0; j < numberOfSteps; j++) {
      movedSpot = movedSpot * exp(pre0 + pre1 * GetOneGaussianByBoxMuller());
    }
    sumPayoff += option.OptionPayOff(movedSpot);
  }
  return (sumPayoff / numberOfPaths);
}

//------------------------------------------------------------------------------
// This function returns the payoff computed using Milstein scheme
//------------------------------------------------------------------------------
double
MonteCarloSimulator::MilsteinWalk(const VanillaOption& option)
{
  double sumPayoff = 0;

  double pre0 = 0.5 * volatility * volatility * deltaT;
  double pre1 = 1 + (riskFreeRate - dividendYield) * deltaT - pre0;
  double pre2 = volatility * sqrt(deltaT);

  for (unsigned long i = 0; i < numberOfPaths; i++) {
    double movedSpot = spot;
    for (unsigned long j = 0; j < numberOfSteps; j++) {
      double gaussian = GetOneGaussianByBoxMuller();
      movedSpot = movedSpot * (pre1 + gaussian * (pre2 + gaussian * pre0));
      /*
      movedSpot = movedSpot * (1 + (riskFreeRate - dividendYield) * deltaT 
                  + volatility * sqrt(deltaT) * gaussian
                  + 0.5 * volatility * volatility * deltaT * (gaussian * gaussian - 1));
      */
    }
    sumPayoff += option.OptionPayOff(movedSpot);
  }

  return (sumPayoff / numberOfPaths);
}
