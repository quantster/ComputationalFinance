
#include <iostream>
#include <cmath>
#include "Simulator.h"
#include "Arrays.h"

//#define ESTIMATE_RHO 1
//#define VANILLA_OPTION_CONTROL 0

using namespace std;

Simulator::Simulator(unsigned long NumberOfPaths_)
  : NumberOfPaths(NumberOfPaths_)
{
  Generator = new RandomRand(NumberOfPaths);
}

Simulator::~Simulator()
{
  delete Generator;
}

void
Simulator::RunSimulation(
          UpAndOutCall& option,
          double Spot,
          double Vol,
          double r,
          double d,
          unsigned long NumberOfDates,
          Statistics& gatherer,
          VariatesType type)
{
  Generator->ResetDimensionality(NumberOfDates);

  // precompute certain values
  double deltaT = option.GetExpiry() / NumberOfDates;

  double pre0 = 1 + (r - d) * deltaT;
  double pre1 = Vol * sqrt(deltaT);

  double discount = exp(-r * option.GetExpiry());

  MJArray variate(NumberOfDates);

  gatherer.Reset();
  switch (type)
  {
    case Normal:
      SimulateWithNormalVariates(option, Spot, pre0, pre1, discount, variate, gatherer);
      break;
    case Antithetic:
      SimulateWithAntitheticVariates(option, Spot, pre0, pre1, discount, variate, gatherer);
      break;
    case Control:
      SimulateWithControlVariates(option, Spot, pre0, pre1, discount, variate, gatherer);
      break;
    default:
      break;
  };

}

void
Simulator::SimulateWithControlVariates(
        UpAndOutCall& option,
        double Spot,
        double pre0,
        double pre1,
        double discount,
        MJArray& variate,
        Statistics& gatherer)
{
  unsigned long NumberOfPilotPaths = 1000;

  double movedSpot;

#ifdef VANILLA_OPTION_CONTROL
    VanillaCall controlOption(option.GetStrike(), option.GetExpiry());
#endif

  // pilot run
  vector<double> x(NumberOfPilotPaths), y(NumberOfPilotPaths);
  double mean_x(0.0), mean_y(0.0);
  for (unsigned long i = 0; i < NumberOfPilotPaths; i++) {

    Generator->GetGaussians(variate);
    movedSpot = Spot;
    option.Reset();

    for (unsigned long j = 0; j < variate.size(); j++) {
      movedSpot = movedSpot * (pre0 + pre1 * variate[j]);

      option.Record(movedSpot);
    }
#ifdef VANILLA_OPTION_CONTROL
    x[i] = controlOption.Payoff(movedSpot);
#else
    x[i] = movedSpot;
#endif
    mean_x += x[i];

    y[i] = discount * option.Payoff(movedSpot);
    mean_y += y[i];
  }


  mean_x = mean_x / NumberOfPilotPaths;
  mean_y = mean_y / NumberOfPilotPaths;

  double numerator(0.0), denominator(0.0);
 
#ifdef ESTIMATE_RHO
  double denominator2(0.0);
#endif

  for (unsigned long i = 0; i < NumberOfPilotPaths; i++) {
    numerator += (x[i] - mean_x) * (y[i] - mean_y);
    denominator += (x[i] - mean_x) * (x[i] - mean_x);
#ifdef ESTIMATE_RHO
    denominator2 += (y[i] - mean_y) * (y[i] - mean_y);
#endif
  }

  // get the estimate for the optimal coefficient
  double b_hat = numerator / denominator;

#ifdef ESTIMATE_RHO
  double rho_hat = numerator / sqrt(denominator * denominator2);
  cout << "b*: " << b_hat << endl;
  cout << "rho: " << rho_hat << endl;
#endif

  for (unsigned long i = 0; i < NumberOfPaths - NumberOfPilotPaths; i++) {
    option.Reset();

    movedSpot = Spot;
    Generator->GetGaussians(variate);

    for (unsigned long j = 0; j < variate.size(); j++) {
      movedSpot = movedSpot * (pre0 + pre1 * variate[j]);

      //if (option.IsKnockedOut(movedSpot)) break;
      option.Record(movedSpot);
    }
    double price = discount * option.Payoff(movedSpot);
    gatherer.DumpOneResult(price - b_hat * (movedSpot - mean_x));
  }
}

void
Simulator::SimulateWithNormalVariates
      (UpAndOutCall& option, 
       double Spot, 
       double pre0, 
       double pre1, 
       double discount, 
       MJArray& variate, 
       Statistics& gatherer)
{
  double movedSpot;
  for (unsigned long i = 0; i < NumberOfPaths; i++) {
    option.Reset();

    movedSpot = Spot;
    Generator->GetGaussians(variate);

    for (unsigned long j = 0; j < variate.size(); j++) {
      movedSpot = movedSpot * (pre0 + pre1 * variate[j]);

      if (option.IsKnockedOut(movedSpot)) break;
    }
    gatherer.DumpOneResult(discount * option.Payoff(movedSpot));
  }
}

void
Simulator::SimulateWithAntitheticVariates
      (UpAndOutCall& option1, 
       double Spot, 
       double pre0, 
       double pre1, 
       double discount, 
       MJArray& variate, 
       Statistics& gatherer)
{
  double movedSpot1, movedSpot2;
  UpAndOutCall option2 = option1;
  for (unsigned long i = 0; i < NumberOfPaths; i++) {
    option1.Reset();
    option2.Reset();

    movedSpot1 = Spot;
    movedSpot2 = Spot;
    Generator->GetGaussians(variate);

    for (unsigned long j = 0; j < variate.size(); j++) {
      movedSpot1 = movedSpot1 * (pre0 + pre1 * variate[j]);
      movedSpot2 = movedSpot1 * (pre0 - pre1 * variate[j]);

      if (option1.IsKnockedOut(movedSpot1) & option2.IsKnockedOut(movedSpot2))
        break;
    }

    gatherer.DumpOneResult(discount * option1.Payoff(movedSpot1), 
              discount * option2.Payoff(movedSpot2));
  }
}

