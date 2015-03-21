// (C) 2010 Nitish Kanabar

#include "Random4.h"
#include "Random1.h"
#include "Normals.h"
#include "Arrays.h"

#include <cmath>

#include <iostream>
using namespace std;

// Fishman Generator
#if UNIFORM_GENERATOR == 1
FishmanGenerator::FishmanGenerator(unsigned long Dimensionality)
  : RandomParkMiller(Dimensionality) {}
#elif UNIFORM_GENERATOR == 2
FishmanGenerator::FishmanGenerator(unsigned long Dimensionality)
  : RandomLEcuyer(Dimensionality) {}
#else
FishmanGenerator::FishmanGenerator(unsigned long Dimensionality)
  : RandomRand(Dimensionality) {}
#endif

void
FishmanGenerator::GetGaussians(MJArray& variates)
{
  for (unsigned long i = 0; i < GetDimensionality(); i++) {
    variates[i] = GetOneGaussian();
  }
}

double
FishmanGenerator::GetOneGaussian()
{
  MJArray u1(1), u2(1), u3(1);
  double x;
  do {
    GetUniforms(u1);
    GetUniforms(u2);
    GetUniforms(u3);
    x = -log(u1[0]);
  }
  while (u2[0] > exp(-0.5 * (x - 1) * (x - 1))); 
  return (u3[0] < 0.5) ? -x : x;
}

// BoxMuller Generator
#if UNIFORM_GENERATOR == 1
BoxMullerGenerator::BoxMullerGenerator(unsigned long Dimensionality)
  : RandomParkMiller(Dimensionality) {}
#elif UNIFORM_GENERATOR == 2
BoxMullerGenerator::BoxMullerGenerator(unsigned long Dimensionality)
  : RandomLEcuyer(Dimensionality) {}
#else
BoxMullerGenerator::BoxMullerGenerator(unsigned long Dimensionality)
  : RandomRand(Dimensionality) {}
#endif

double
BoxMullerGenerator::GetOneGaussian()
{
  MJArray u1(1);
  MJArray u2(1);

  double x, y;
  do {
    GetUniforms(u1);
    GetUniforms(u2);

    u1[0] = 2 * u1[0] - 1;
    u2[0] = 2 * u2[0] - 1;
    x = u1[0]*u1[0] + u2[0]*u2[0];
    y = sqrt(-2 * log(x) / x);
  }
  while (x > 1);

  return u1[0] * y;
}

void
BoxMullerGenerator::GetGaussians(MJArray& variates)
{
  for (unsigned long i = 0; i < GetDimensionality(); i++) {
    variates[i] = GetOneGaussian();
  }
}


// InverseTransform Generator
#if UNIFORM_GENERATOR == 1
InverseTransformGenerator::InverseTransformGenerator(unsigned long Dimensionality)
  : RandomParkMiller(Dimensionality) {}
#elif UNIFORM_GENERATOR == 2
InverseTransformGenerator::InverseTransformGenerator(unsigned long Dimensionality)
  : RandomLEcuyer(Dimensionality) {}
#else
InverseTransformGenerator::InverseTransformGenerator(unsigned long Dimensionality)
  : RandomRand(Dimensionality) {}
#endif

double
InverseTransformGenerator::GetOneGaussian()
{
  MJArray uniforms(1);
  GetUniforms(uniforms);
  return InverseCumulativeNormal(uniforms[0]);
}

void
InverseTransformGenerator::GetGaussians(MJArray& variates)
{
  for (unsigned long i = 0; i < GetDimensionality(); i++) {
    variates[i] = GetOneGaussian();
  }
}

