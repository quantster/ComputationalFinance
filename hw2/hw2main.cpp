#include <iostream>
#include <sys/time.h>

#include "Simulator.h"
#include "Option.h"

using namespace std;

//
// Function get_time as implemented by Dr. Paul Feehan
//
double get_time()
{
    struct timeval t;
    gettimeofday(&t, NULL);
    double d = t.tv_sec + (double) t.tv_usec/1000000;
    return d;
}

int main(int argc, char* argv[])
{
	double Spot(100);
	double r(0.05);
	double d(0.02);
	double Vol(0.3);

	double Strike(110);
  double Expiry(1);

	double Barrier(120);
  UpAndOutCall option(Strike, Barrier, Expiry);

	unsigned long NumberOfPaths(10000);
	unsigned long NumberOfDates(252);

  Statistics gatherer;
  Simulator engine(NumberOfPaths);

  double start_t, end_t;
  vector<double> results;
  start_t = get_time();
  double price = option.Price(Spot, Vol, r, d); 
  end_t = get_time();
  cout << "Closed-form barrier option price = " << price
       << ", run time = " << end_t - start_t << " seconds." << endl;

  start_t = get_time();
  engine.RunSimulation(option, Spot, Vol, r, d, NumberOfDates, gatherer, Simulator::Normal);
  results = gatherer.GetResultsSoFar();
  end_t = get_time();
  cout << "Monte Carlo barrier option price = " << results[0]
       << ", std error = " << results[1]
       << ", run time = " << end_t - start_t << " seconds." << endl;

  start_t = get_time();
  engine.RunSimulation(option, Spot, Vol, r, d, NumberOfDates, gatherer, Simulator::Antithetic);
  results = gatherer.GetResultsSoFar();
  end_t = get_time();
  cout << "Antithetic variates barrier option price = " << results[0]
       << ", std error = " << results[1]
       << ", run time = " << end_t - start_t << " seconds." <<  endl;

  start_t = get_time();
  engine.RunSimulation(option, Spot, Vol, r, d, NumberOfDates, gatherer, Simulator::Control);
  results = gatherer.GetResultsSoFar();
  end_t = get_time();
  cout << "Control Variates barrier option price = " << results[0]
       << ", std error = " << results[1]
       << ", run time = " << end_t - start_t << " seconds." << endl;


  return 0;
}
