#include "LEcuyer.h"

#include <iostream>

using namespace std;

MultipleRecursiveGenerator::MultipleRecursiveGenerator(
      unsigned long _m,
      unsigned long _a1,
      unsigned long _a2,
      unsigned long _a3,
      long _seed1,
      long _seed2,
      long _seed3)
  : m(_m), a1(_a1), a2(_a2), a3(_a3), seed1(_seed1), seed2(_seed2), seed3(_seed3)
{
}

long
MultipleRecursiveGenerator::GetOneRandomInteger()
{
  long random = (a1 * seed3 + a2 * seed2 + a3 * seed1) % m;
  seed1 = seed2;
  seed2 = seed3;
  seed3 = random;
  return random;
}

RandomLEcuyer::RandomLEcuyer(unsigned long Dimensionality)
  : RandomRand(Dimensionality)
{
  mrg1 = new MultipleRecursiveGenerator(2147483647, 0, 63308, -183326, 1, 2, 3);
  mrg2 = new MultipleRecursiveGenerator(2145483479, 86098, 0, -539608, 5, 7, 11);
}

RandomLEcuyer::~RandomLEcuyer()
{
  delete mrg1;
  delete mrg2;
}


double
RandomLEcuyer::GetOneUniform()
{
  long random1 = mrg1->GetOneRandomInteger();
  long random2 = mrg2->GetOneRandomInteger();
  
  long random = (random1 - random2) % mrg1->M();

  return random / double(mrg1->M());
}

void
RandomLEcuyer::GetUniforms(MJArray& variates)
{
  for (unsigned long j=0; j < GetDimensionality(); j++)
    variates[j] = GetOneUniform();
}
