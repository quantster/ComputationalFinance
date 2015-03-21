//
//
//		RandomMain3.cpp
//
//       
//

//#define PRINT_GAUSSIANS

#include <iostream>

#include <SimpleMC8.h>
#include <ParkMiller.h>
#include <Vanilla3.h>
#include <MCStatistics.h>
#include <BlackScholesFormulas.h>
#include "LEcuyer.h"
#include "Random4.h"

using namespace std;

int main()
{
	double Expiry(1);
	double Strike(110); 
	double Spot(100); 
	double Vol(0.3); 
	double r(0.05); 
	unsigned long NumberOfPaths(10000);


  PayOffCall thePayOff(Strike);

  VanillaOption theOption(thePayOff, Expiry);

  ParametersConstant VolParam(Vol);
  ParametersConstant rParam(r);

  vector<vector<double> > results;

  // closed-form solution
  {
    cout << "Closed-form option price = "
        << BlackScholesCall(Spot, Strike, r, 0.0, Vol, Expiry) << endl;
  }

  // Pricing using Park-Miller Uniform generator
  {
    StatisticsMean gatherer;
    RandomParkMiller generator(1);
//    AntiThetic GenTwo(generator);

	  SimpleMonteCarlo6(theOption,
        Spot, 
        VolParam,
        rParam,
        NumberOfPaths,
        gatherer,
        generator);

    results = gatherer.GetResultsSoFar();

    cout << "MC option price with Park-Miller uniform generator = "
         << results[0][0] << endl;
  }

  // pricing using LEcuyer uniform generator
  {
    StatisticsMean gatherer;

    RandomLEcuyer generator(1);

	  SimpleMonteCarlo6(theOption,
        Spot, 
        VolParam,
        rParam,
        NumberOfPaths,
        gatherer,
        generator);

    results = gatherer.GetResultsSoFar();

    cout << "MC option price with L'Ecuyer uniform generator = "
         << results[0][0] << endl;
  }

  // pricing using inverse distribution normal generation
  {
    StatisticsMean gatherer;
    InverseTransformGenerator generator(1);

	  SimpleMonteCarlo6(theOption,
        Spot, 
        VolParam,
        rParam,
        NumberOfPaths,
        gatherer,
        generator);

    results = gatherer.GetResultsSoFar();
    cout << "MC option price with inverse distribution normal generator = "
         << results[0][0] << endl;
  }

  // Box-Muller normal generator
  {
    StatisticsMean gatherer;
    BoxMullerGenerator generator(1);

	  SimpleMonteCarlo6(theOption,
        Spot, 
        VolParam,
        rParam,
        NumberOfPaths,
        gatherer,
        generator);

    results = gatherer.GetResultsSoFar();
    cout << "MC option price with Box-Muller normal generator = "
         << results[0][0] << endl;

  }

  // Fishman normal generator
  {
    StatisticsMean gatherer;
    FishmanGenerator generator(1);

	  SimpleMonteCarlo6(theOption,
        Spot, 
        VolParam,
        rParam,
        NumberOfPaths,
        gatherer,
        generator);

    results = gatherer.GetResultsSoFar();
    cout << "MC option price with Fishman normal generator = "
         << results[0][0] << endl;

  }

  // Question 2
#ifdef PRINT_GAUSSIANS
#if UNIFORM_GENERATOR == 0
  cout << "Using RandomRand Uniforms" << endl;
#elif UNIFORM_GENERATOR == 1
  cout << "Using ParkMiller Uniforms" << endl;
#else
  cout << "Using LEcuyer Uniforms" << endl;
#endif
  unsigned long numberOfSamples(1000);

  //cout << "InverseTransform_Gaussians" << endl;
  InverseTransformGenerator inverseTransform(1);
  BoxMullerGenerator boxMuller(1);
  FishmanGenerator fishman(1);

  cout << "Inverse,BoxMuller,Fishman" << endl;
  for (unsigned long i = 0; i < numberOfSamples; i++) {
    cout << inverseTransform.GetOneGaussian() << ", "
        << boxMuller.GetOneGaussian() << ", "
        << fishman.GetOneGaussian() << endl;
  }

#endif

  return 0;
}

/*
 *
 * Copyright (c) 2002
 * Mark Joshi
 *
 * Permission to use, copy, modify, distribute and sell this
 * software for any purpose is hereby
 * granted without fee, provided that the above copyright notice
 * appear in all copies and that both that copyright notice and
 * this permission notice appear in supporting documentation.
 * Mark Joshi makes no representations about the
 * suitability of this software for any purpose. It is provided
 * "as is" without express or implied warranty.
*/

