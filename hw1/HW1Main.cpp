////////////////////////////////////////////////////////////////////////////////
//
//  Author: Nitish Kanabar
//  Course: Math 623, Spring 2010
//  Assignment: HW 1
//  Compiler Used: 686-apple-darwin10-g++-4.2.1 (GCC) 4.2.1 (Apple Inc. build 5646) (dot 1)
//
////////////////////////////////////////////////////////////////////////////////

#define DEBUG 1

#include <iostream>
#include <cmath>
#include <BlackScholesFormulas.h>
#include "MonteCarloSimulator.h"

using namespace std;


int main(int argc, char* argv[])
{
  // default option parameters
  double spot(100);
  double riskFreeRate(0.05);
  double volatility(0.3);
  double dividendYield(0.02);
  double strike(110);
  double expiry(1);

  enum OptionType{call, put};

  // defaults for the simulation
  unsigned long numSteps(252);
  unsigned long numPaths(10000);

  // MonteCarloSimulator(Spot, Rate, Volatility, DividendYield, NumPaths, NumSteps)
  MonteCarloSimulator pricer(spot, riskFreeRate, volatility, dividendYield, numPaths, numSteps);

  // defaults for the assignment
  OptionType optionType(call);

  PayOffCall callPayoff(strike);
  VanillaOption vanillaCall(callPayoff, expiry);

  PayOffPut putPayoff(strike);
  VanillaOption vanillaPut(putPayoff, expiry);

  // closed-form pricing
  cout << "Option price using closed-form formula = ";
  if (optionType == call)
    cout << BlackScholesCall(spot, strike, riskFreeRate, dividendYield, volatility, expiry);
  else
    cout << BlackScholesPut(spot, strike, riskFreeRate, dividendYield, volatility, expiry);
  cout << endl;

  // single step exact SDE solution
  cout << "Option price using single-step exact SDE solution = ";
  if (optionType == call)
    cout << pricer.ExpectedPayoff(vanillaCall, MonteCarloSimulator::SingleStep);
  else
    cout << pricer.ExpectedPayoff(vanillaPut, MonteCarloSimulator::SingleStep);
  cout << endl;

  // Euler spot
  cout << "Option price using Euler numerical solution of SDE for spot = ";
  if (optionType == call)
    cout << pricer.ExpectedPayoff(vanillaCall, MonteCarloSimulator::EulerSpot);
  else
    cout << pricer.ExpectedPayoff(vanillaPut, MonteCarloSimulator::EulerSpot);
  cout << endl;

  // Euler log spot
  cout << "Option price using Euler numerical solution of SDE for log spot = ";
  if (optionType == call)
    cout << pricer.ExpectedPayoff(vanillaCall, MonteCarloSimulator::EulerLogSpot);
  else
    cout << pricer.ExpectedPayoff(vanillaPut, MonteCarloSimulator::EulerLogSpot);
  cout << endl;

  // Milstein
  cout << "Option price using Milstein numerical solution of SDE for spot = ";
  if (optionType == call)
    cout << pricer.ExpectedPayoff(vanillaCall, MonteCarloSimulator::Milstein);
  else
    cout << pricer.ExpectedPayoff(vanillaPut, MonteCarloSimulator::Milstein);
  cout << endl;

  return 0;
}
