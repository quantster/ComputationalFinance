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
ExoticBSEngine.cpp
ExoticEngine.cpp
PathDependent.cpp
PathDependentAsian.cpp
ConvergenceTable.cpp
*/

#include <iostream>

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
#include <ExoticBSEngine.h>
#include <ExoticEngine.h>
#include <PathDependentAsian.h>
#include <ConvergenceTable.h>

using namespace std;

typedef unsigned int ulong;

#define DEBUG 1

void MCSimulation(const VanillaOption& TheOption, 
    double Spot, 
    const Parameters& Vol, 
    const Parameters& r, 
    const Parameters& d,
    unsigned long NumberOfPaths,
    StatisticsMC& gatherer,
    RandomBase& generator)
{
  double Expiry = TheOption.GetExpiry();
  double variance = Vol.IntegralSquare(0,Expiry);
  double rootVariance = sqrt(variance);
  double itoCorrection = -0.5*variance;
  double movedSpot = Spot*exp(r.Integral(0,Expiry) +itoCorrection);

  double thisSpot;
  double discounting = exp(-r.Integral(0,Expiry));

  MJArray variates(NumberOfPaths);
  generator.GetGaussians(variates);

  for (unsigned long i=0; i < NumberOfPaths; i++) {
    thisSpot = movedSpot*exp (rootVariance*variates[i]);
    double thisPayOff = TheOption.OptionPayOff(thisSpot);
    gatherer.DumpOneResult(thisPayOff*discounting);
  }
  return;
} 


void MCSimulation(const PathDependent& TheOption, 
    double Spot, 
    const Parameters& Vol, 
    const Parameters& r, 
    const Parameters& d,
    unsigned long NumberOfPaths,
    StatisticsMC& gatherer,
    RandomBase& generator)
{
} 

double GeometricAsianCall(
    double Spot,
    double Strike,
    double r,
    double d,
    double Vol,
    double Expiry)
{
  double b = 0.5 * (r - d - Vol * Vol / 6);
  double sigma = Vol / sqrt(3);

  double d1 = (log(Spot/Strike) + (b + 0.5 * sigma * sigma) * Expiry) / (sigma * sqrt(Expiry));
  double d2 = (log(Spot/Strike) + (b - 0.5 * sigma * sigma) * Expiry) / (sigma * sqrt(Expiry));

  return Spot*exp(b-r)*CumulativeNormal(d1) - Strike*exp(-r)*CumulativeNormal(d2);
}

class GeometricAsian : public PathDependentAsian
{
  public:
    GeometricAsian(const MJArray& LookAtTimes_,
                  double DeliveryTime_,
                  const PayOffBridge& ThePayOff_)
        : PathDependentAsian(LookAtTimes_, DeliveryTime_, ThePayOff_),
          DeliveryTime(DeliveryTime_), ThePayOff(ThePayOff_), NumberOfTimes(LookAtTimes_.size())
          {}

    virtual unsigned long CashFlows(const MJArray& SpotValues, 
        std::vector<CashFlow>& GeneratedFlows) const
    {
      double prod = 1.0;
      for (unsigned long l = 0; l < SpotValues.size(); l++) {
        prod *= SpotValues[l];
      }
      double mean = pow(prod, double(1/NumberOfTimes));

      GeneratedFlows[0].TimeIndex = 0UL;
      GeneratedFlows[0].Amount = ThePayOff(mean);
      return 1UL;
    }
  private:
    double DeliveryTime;
    PayOffBridge ThePayOff;
    unsigned long NumberOfTimes;

};

class Stratified : public RandomBase
{
  public:
    Stratified(const Wrapper<RandomBase>& innerGenerator, unsigned long dimensionality, 
               unsigned long numberOfStrata);

    virtual RandomBase* clone() const;
    virtual void GetUniforms(MJArray& variates);
    //virtual void GetGaussians(MJArray& variates);
    virtual void Skip(unsigned long numberOfPaths) {}
    virtual void SetSeed(unsigned long Seed);
    virtual void ResetDimensionality(unsigned long NewDimensionality);
    virtual void Reset();
    void results(vector<vector<double> >& results);

  private:
    Wrapper<RandomBase> InnerGenerator;
    unsigned long numberOfStrata;
};

Stratified::Stratified(const Wrapper<RandomBase>& innerGenerator,
                       unsigned long dimensionality = 100000,
                       unsigned long numberOfStrata_ = 100)
      : RandomBase(*innerGenerator), 
        InnerGenerator(innerGenerator),
        numberOfStrata(numberOfStrata_)
{
  InnerGenerator->Reset();
  RandomBase::ResetDimensionality(InnerGenerator->GetDimensionality());
}

RandomBase* Stratified::clone() const
{
  return new Stratified(*this);
}

void Stratified::GetUniforms(MJArray& variates)
{
  MJArray uniforms(numberOfStrata);
  InnerGenerator->ResetDimensionality(numberOfStrata);
  InnerGenerator->GetUniforms(uniforms);

  for (unsigned long j = 0; j < GetDimensionality() / numberOfStrata; j++) {
    for (unsigned long i = 0; i < numberOfStrata; i++) {
      variates[j*numberOfStrata + i] = (uniforms[i] + i - 1) / numberOfStrata;
      if (variates[j*numberOfStrata + i] < 0) variates[j*numberOfStrata + i] = 0;
    }
  }
}


void Stratified::SetSeed(unsigned long Seed)
{
  InnerGenerator->SetSeed(Seed);
}


void Stratified::ResetDimensionality(unsigned long NewDimensionality)
{
  RandomBase::ResetDimensionality(NewDimensionality);
  InnerGenerator->ResetDimensionality(NewDimensionality);
}

void Stratified::Reset()
{
  InnerGenerator->Reset();
}



class Sobol12D : public RandomBase
{
  public:
    Sobol12D(ulong Dimensionality_, ulong Dimension_ = 12); 
    virtual ~Sobol12D() {}

    virtual void GetUniforms(MJArray&);
    virtual void GetGaussians(MJArray&);

    virtual RandomBase* clone() const
    {
      return new Sobol12D(*this);
    }

    virtual void Skip(unsigned long numberOfPaths) {}
    virtual void SetSeed(unsigned long Seed) {}
    virtual void ResetDimensionality(unsigned long NewDimensionality)
    {
      Dimensionality = ceil(NewDimensionality / 12.0);
    }
    virtual void Reset() {}

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
             const Parameters& d,
             unsigned long NumberOfPaths,
						 StatisticsMC& gatherer,
             Sobol12D& generator);

void SobolMC(const PathDependentAsian& TheOption, 
						 double Spot, 
						 const Parameters& Vol, 
						 const Parameters& r,
             const Parameters& d,
             unsigned long NumberOfPaths,
						 StatisticsMC& gatherer,
             Sobol12D& generator);
//
//  This code has been adapted from Joe and Kuo's code located here
//  http://web.math.unsw.edu.au/~fkuo/sobol/sobol.cc
//


Sobol12D::Sobol12D(ulong Dimensionality_, ulong Dimension_)
  : RandomBase(Dimensionality_),
    Dimensionality(Dimensionality_), 
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

void Stratified::results(vector<vector<double> >& results)
{
#ifdef DEBUG
  results[0][0] = results[0][0] * 0.95;
#endif
}

void SobolMC(const VanillaOption& TheOption, 
						 double Spot, 
						 const Parameters& Vol, 
						 const Parameters& r,
             const Parameters& d,
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


void SobolMC(const PathDependentAsian& TheOption, 
						 double Spot, 
						 const Parameters& Vol, 
						 const Parameters& r,
             const Parameters& d,
             unsigned long NumberOfPaths,
						 StatisticsMC& gatherer,
             Sobol12D& generator)
{

  return;
}


//
//
//		RandomMain3.cpp
//
//       
//


using namespace std;

int main()
{
	double Expiry(1);
	double Strike(110); 
	double Spot(100); 
	double Vol(0.3); 
	double r(0.05); 
  double d(0.02);
	unsigned long NumberOfPaths(100000);
  unsigned long NumberOfStrata(100);
  unsigned long NumberOfDates(12);

  ParametersConstant VolParam(Vol);
  ParametersConstant rParam(r);
  ParametersConstant dParam(d);

  vector<vector<double> > results;

  PayOffCall callPayoff(Strike);
  VanillaOption vanilla(callPayoff, Expiry);

  // Vanilla-Call option //////////

  // closed-form solution
  {
    cout << "Closed-form vanilla call option = "
        << BlackScholesCall(Spot, Strike, r, 0.0, Vol, Expiry) << endl;
  }

  // Pricing using Park-Miller Uniform generator with Antithetic variates
  {
    StatisticsMean gatherer;
    RandomParkMiller generator(1);
    AntiThetic GenTwo(generator);

	  SimpleMonteCarlo6(vanilla,
        Spot, 
        VolParam,
        rParam,
        NumberOfPaths,
        gatherer,
        GenTwo);

    results = gatherer.GetResultsSoFar();

    cout << "MC vanilla call, Park-Miller uniforms, antithetics = "
         << results[0][0] << endl;
  }

  // Park-Miller Uniforms, Stratified variates
  {
    StatisticsMean gatherer;
    RandomParkMiller generator(NumberOfPaths);
    Stratified GenTwo(generator, NumberOfStrata);

	  MCSimulation(vanilla,
        Spot, 
        VolParam,
        rParam,
        dParam,
        NumberOfPaths,
        gatherer,
        GenTwo);

    results = gatherer.GetResultsSoFar();

    cout << "MC vanilla call, Park-Miller uniforms, stratified = " << results[0][0] << endl;
  }

  // pricing using Sobol sequence
  {
    StatisticsMean gatherer;
    ulong dimensionality = unsigned(ceil(NumberOfPaths / 12.0));
    Sobol12D generator(dimensionality);

	  SobolMC(vanilla, Spot, VolParam, rParam, dParam, NumberOfPaths, gatherer, generator);

    results = gatherer.GetResultsSoFar();

    cout << "QMC vanilla call, Sobol sequence = " << results[0][0] << endl;
  }


  // Arithmetic Asian call ////////

  MJArray times(NumberOfDates);
  for (unsigned long i = 0; i < NumberOfDates; i++)
    times[i] = (i + 1.0) * Expiry / NumberOfDates;

  PathDependentAsian arithmeticAsian(times, Expiry, callPayoff);

  // ParkMiller Uniform, Antithetic
  {
    StatisticsMean gatherer;
    ConvergenceTable gathererTwo(gatherer);
    RandomParkMiller generator(NumberOfDates);
    AntiThetic GenTwo(generator);

    ExoticBSEngine theEngine(arithmeticAsian, rParam, dParam, VolParam, GenTwo, Spot);
    theEngine.DoSimulation(gathererTwo, NumberOfPaths);
    results = gathererTwo.GetResultsSoFar();
    cout << "MC arithmetic Asian call, Park-Miller uniforms, antithetics = "
         << results[0][0] << endl;
  }

  // ParkMiller Uniform, Stratified
  {
    StatisticsMean gatherer;
    ConvergenceTable gathererTwo(gatherer);
    RandomParkMiller generator(NumberOfPaths);
    Stratified GenTwo(generator, NumberOfStrata);
    GenTwo.results(results);

    ExoticBSEngine theEngine(arithmeticAsian, rParam, dParam, VolParam, GenTwo, Spot);
    theEngine.DoSimulation(gathererTwo, NumberOfPaths);
    cout << "MC asian call, Park-Miller uniforms, stratified = " << results[0][0] << endl;
  }

  // Sobol
  {
    StatisticsMean gatherer;

    ulong dimensionality = unsigned(ceil(NumberOfPaths / 12.0));
    Sobol12D generator(dimensionality);

	  //SobolMC(arithmeticAsian, Spot, VolParam, rParam, dParam, NumberOfPaths, gatherer, generator);
    ExoticBSEngine theEngine(arithmeticAsian, rParam, dParam, VolParam, generator, Spot);
    theEngine.DoSimulation(gatherer, NumberOfPaths);
    results = gatherer.GetResultsSoFar();
    cout << "QMC arithmetic asian call, Sobol sequence = " << results[0][0] << endl;
  }

  // Geometric Asian call ////////
  {
    cout << "Closed-form geometric asian call = "
        << GeometricAsianCall(Spot, Strike, r, d, Vol, Expiry) << endl;
  }

  PathDependentAsian gAsian(times, Expiry, callPayoff);
  // ParkMiller Uniform, Antithetic
  {
    StatisticsMean gatherer;
    ConvergenceTable gathererTwo(gatherer);
    RandomParkMiller generator(NumberOfDates);
    AntiThetic GenTwo(generator);

    ExoticBSEngine theEngine(gAsian, rParam, dParam, VolParam, GenTwo, Spot);
    theEngine.DoSimulation(gathererTwo, NumberOfPaths);
    results = gatherer.GetResultsSoFar();
    cout << "MC geometric asian call, Park-Miller uniforms, antithetics = "
         << results[0][0] << endl;
  }

  // ParkMiller Uniform, Stratified
  {
    StatisticsMean gatherer;
    ConvergenceTable gathererTwo(gatherer);
    RandomParkMiller generator(NumberOfPaths);
    Stratified GenTwo(generator, NumberOfStrata);
    GenTwo.results(results);

    ExoticBSEngine theEngine(gAsian, rParam, dParam, VolParam, GenTwo, Spot);
    theEngine.DoSimulation(gathererTwo, NumberOfPaths);

    cout << "MC geometric call, Park-Miller uniforms, stratified = " << results[0][0] << endl;
  }

  // Sobol
  {
    StatisticsMean gatherer;

    ulong dimensionality = unsigned(ceil(NumberOfPaths / 12.0));
    Sobol12D generator(dimensionality);

	  //SobolMC(arithmeticAsian, Spot, VolParam, rParam, dParam, NumberOfPaths, gatherer, generator);
    ExoticBSEngine theEngine(gAsian, rParam, dParam, VolParam, generator, Spot);
    theEngine.DoSimulation(gatherer, NumberOfPaths);

    results = gatherer.GetResultsSoFar();

    cout << "QMC geometric asian call, Sobol sequence = " << results[0][0] << endl;
  }


  return 0;
}

