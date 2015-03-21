//
//  This code has been adapted from Joe and Kuo's code located here
//  http://web.math.unsw.edu.au/~fkuo/sobol/sobol.cc
//

/*
The following files need to be compiled
Arrays.cpp
Vanilla3.cpp
Parameters.cpp
Random2.cpp
MCStatistics.cpp
Normals.cpp
SimpleMC8.cpp
ParkMiller.cpp
BlackScholesFormulas.cpp
Antithetic.cpp
*/

#include <vector>
#include <cmath>
#include <Arrays.h>
#include <Vanilla3.h>
#include <Parameters.h>
#include <Random2.h>
#include <MCStatistics.h>
#include <Normals.h>
#include <iostream>
#include <SimpleMC8.h>
#include <ParkMiller.h>
#include <BlackScholesFormulas.h>
#include <Antithetic.h>

using namespace std;

typedef unsigned int ulong;

class Sobol12D
{
  public:
    Sobol12D(ulong Dimensionality_, ulong Dimension_ = 12); 
    virtual ~Sobol12D() {}

    void GetUniforms(MJArray&);
    void GetGaussians(MJArray&);


  private:
    ulong Dimensionality;
    ulong Dimension;

    vector<vector<ulong> > m;
    vector<ulong> a;
};


void SobolMC(const VanillaOption& TheOption, 
						 double Spot, 
						 const Parameters& Vol, 
						 const Parameters& r,
             unsigned long NumberOfPaths,
						 StatisticsMC& gatherer,
             Sobol12D& generator);
//
//  This code has been adapted from Joe and Kuo's code located here
//  http://web.math.unsw.edu.au/~fkuo/sobol/sobol.cc
//


Sobol12D::Sobol12D(ulong Dimensionality_, ulong Dimension_)
  : Dimensionality(Dimensionality_), 
    Dimension(Dimension_),
    m(vector<vector<ulong> >(11)), 
    a(vector<ulong>(11))
{
  // dimension 2
  a[0] = 0;
  m[0] = vector<ulong>(1);
  m[0][0] = 1;

  // dimension 3
  a[1] = 1;
  m[1] = vector<ulong>(2);
  m[1][0] = 1; m[1][1] = 1;

  // dimension 4
  a[2] = 1;
  m[2] = vector<ulong>(3);
  m[2][0] = 1; m[2][1] = 1; m[2][2] = 1;

  // dimension 5
  a[3] = 2;
  m[3] = vector<ulong>(3);
  m[3][0] = 1; m[3][1] = 3; m[3][2] = 1;

  // dimension 6
  a[4] = 1;
  m[4] = vector<ulong>(4);
  m[4][0] = 1; m[4][1] = 1; m[4][2] = 7; m[4][3] = 13;

  // dimension 7
  a[5] = 4;
  m[5] = vector<ulong>(4);
  m[5][0] = 1; m[5][1] = 1; m[5][2] = 3; m[5][3] = 7;

  // dimension 8
  a[6] = 2;
  m[6] = vector<ulong>(5);
  m[6][0] = 1; m[6][1] = 3; m[6][2] = 1; m[6][3] = 7; m[6][4] = 21;

  // dimension 9
  a[7] = 13;
  m[7] = vector<ulong>(5);
  m[7][0] = 1; m[7][1] = 3; m[7][2] = 1; m[7][3] = 3; m[7][4] = 9;

  // dimension 10
  a[8] = 7;
  m[8] = vector<ulong>(5);
  m[8][0] = 1; m[8][1] = 1; m[8][2] = 5; m[8][3] = 9; m[8][4] = 13;

  // dimension 11
  a[9] = 14;
  m[9] = vector<ulong>(5);
  m[9][0] = 1; m[9][1] = 1; m[9][2] = 3; m[9][3] = 9; m[9][4] = 13;

  // dimension 12
  a[10] = 11;
  m[10] = vector<ulong>(5);
  m[10][0] = 1; m[10][1] = 1; m[10][2] = 5; m[10][3] = 3; m[10][4] = 7;
}

void
Sobol12D::GetUniforms(MJArray& v)
{

  ulong NumberOfBits = unsigned(ceil(log(double(Dimensionality))/log(2.0)));

  // index from the right to the first zero bit of i
  ulong C[Dimensionality];
  C[0] = 1;
  for (ulong i = 1; i <= Dimensionality - 1; i++) {
    C[i] = 1;
    ulong value = i;
    while (value & 1) {
      value >>= 1;
      C[i]++;
    }
  }

  vector<vector<double> > variates(Dimensionality + 1);
  for (ulong i = 0; i <= Dimensionality - 1; i++) {
    variates[i] = vector<double>(Dimension);
    for (ulong j = 0; j <= Dimension - 1; j++) 
      variates[0][j] = 0;
  }


  // compute the first dimension /////////////////////

  // compute direction numbers V[1] to V[L], scaled by pow(2, 32)
  ulong V[NumberOfBits + 1];
  for (ulong i = 1; i <= NumberOfBits; i++) {
    V[i] = 1 << (32 - i);
  }


  // Evalulate X[0] to X[N-1], scaled by pow(2,32)
  ulong X[Dimensionality];
  X[0] = 0;
  for (ulong i = 1; i <= Dimensionality - 1; i++) {
    X[i] = X[i - 1] ^ V[C[i - 1]];
    variates[i][0] = double(X[i])/pow(2.0,32); 
  }


  // compute the remaining dimensions ////////////////
  for (ulong j = 1; j <= Dimension - 1; j++) {
    ulong s = m[j-1].size();

    ulong V[NumberOfBits + 1];
    if (NumberOfBits <= s) {
      for (ulong i = 1; i <= NumberOfBits; i++) {
        V[i] = m[j-1][i-1] << (32-i); 
      }
    }
    else {
      for (ulong i = 1; i <= s; i++)
        V[i] = m[j-1][i-1] << (32-i); 

      for (ulong i = s + 1; i <= NumberOfBits; i++) {
        V[i] = V[i - s] ^ (V[i - s] >> s); 

        for (ulong k = 1; k <= s - 1; k++) 
          V[i] ^= (((a[j-1] >> (s-1-k)) & 1) * V[i-k]); 
      }
    }

    // Evalulate X[0] to X[N-1], scaled by pow(2,32)
    ulong X[Dimensionality];
    X[0] = 0;
    for (ulong i = 1; i <= Dimensionality - 1; i++) {
      X[i] =  X[i-1] ^ V[C[i-1]];
      variates[i][j] = double(X[i])/pow(2.0,32);
      //          ^ j for dimension (j+1)
    }
  }

  /*
  for (ulong i = 1; i < Dimensionality; i++) {
    for (ulong j = 0; j < Dimension; j++) {
      v[(i-1) * Dimension + j] = variates[i][j];
    }
  }
  */
  for (ulong i = 0; i < Dimensionality; i++) {
    for (ulong j = 0; j < Dimension; j++) {
      v[i * Dimension + j] = variates[i][j];
    }
  }
}

void
Sobol12D::GetGaussians(MJArray& variates)
{
  GetUniforms(variates);
  for (ulong i = 0; i < Dimensionality * Dimension; i++)
    variates[i] = InverseCumulativeNormal(variates[i]);
}

using namespace std;

void SobolMC(const VanillaOption& TheOption, 
						 double Spot, 
						 const Parameters& Vol, 
						 const Parameters& r,
             unsigned long NumberOfPaths,
						 StatisticsMC& gatherer,
             Sobol12D& generator)
{
  double Expiry = TheOption.GetExpiry();
  double variance = Vol.IntegralSquare(0,Expiry);
  double rootVariance = sqrt(variance);
  double itoCorrection = -0.5*variance;
  double movedSpot = Spot*exp(r.Integral(0,Expiry) +itoCorrection);

  double thisSpot;
  double discounting = exp(-r.Integral(0,Expiry));

  MJArray variates(NumberOfPaths + 1);
  generator.GetGaussians(variates);

  for (unsigned long i = 0; i < NumberOfPaths; i++) {
    thisSpot = movedSpot * exp(rootVariance * variates[i]);
    double thisPayOff = TheOption.OptionPayOff(thisSpot);
    gatherer.DumpOneResult(thisPayOff * discounting);
  }

  return;
}


//
//
//		RandomMain3.cpp
//
//       
//

//#define PRINT_UNIFORMS


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
    cout << "Closed-form vanilla call option price = "
        << BlackScholesCall(Spot, Strike, r, 0.0, Vol, Expiry) << endl;
  }

  // Pricing using Park-Miller Uniform generator
  {
    StatisticsMean gatherer;
    RandomParkMiller generator(1);

	  SimpleMonteCarlo6(theOption,
        Spot, 
        VolParam,
        rParam,
        NumberOfPaths,
        gatherer,
        generator);

    results = gatherer.GetResultsSoFar();

    cout << "MC vanilla call option price with Park-Miller uniforms = "
         << results[0][0] << endl;
  }

  // Pricing using Park-Miller Uniform generator with Antithetic variates
  {
    StatisticsMean gatherer;
    RandomParkMiller generator(1);
    AntiThetic GenTwo(generator);

	  SimpleMonteCarlo6(theOption,
        Spot, 
        VolParam,
        rParam,
        NumberOfPaths,
        gatherer,
        GenTwo);

    results = gatherer.GetResultsSoFar();

    cout << "MC vanilla call option price with Park-Miller uniforms and antithetics = "
         << results[0][0] << endl;
  }

  // pricing using Sobol sequence
  {
    StatisticsMean gatherer;

    ulong dimensionality = unsigned(ceil(NumberOfPaths / 12.0));
    Sobol12D generator(dimensionality);

	  SobolMC(theOption,
        Spot, 
        VolParam,
        rParam,
        NumberOfPaths,
        gatherer,
        generator);

    results = gatherer.GetResultsSoFar();

    cout << "QMC vanilla call price with Sobol sequence = "
         << results[0][0] << endl;

  }

  // Question 2
#ifdef PRINT_UNIFORMS
  // 1024 Park-Miller pairs
  RandomParkMiller parkmiller(2);
  Sobol12D sobol(1024, 2);

  MJArray sobol_variates(1024 * 2);
  MJArray pm_variates(2);

  sobol.GetUniforms(sobol_variates);
  for (unsigned int i = 0; i <= 1024; i++) {
    parkmiller.GetUniforms(pm_variates);
    cout << sobol_variates[i] << "," << sobol_variates[++i] << ","
         << pm_variates[0] << "," << pm_variates[1] << endl;
  }

#endif

  return 0;
}

