//  Nitish Kanabar

#ifndef LECUYER_H
#define LECUYER_H

#include "Random3.h"
#include "Arrays.h"

class MultipleRecursiveGenerator
{
  public:
    MultipleRecursiveGenerator(
          unsigned long _m,
          unsigned long _a1,
          unsigned long _a2,
          unsigned long _a3,
          long _seed1,
          long _seed2,
          long _seed3);

    long GetOneRandomInteger();
    inline unsigned long M() const { return m; }

  private:
    unsigned long m;
    unsigned long a1;
    unsigned long a2;
    unsigned long a3;
    long seed1;
    long seed2;
    long seed3;
};

class RandomLEcuyer : public RandomRand
{
  public:
    RandomLEcuyer(unsigned long Dimensionality);
    virtual ~RandomLEcuyer();

    virtual void GetUniforms(MJArray& variates);
    virtual double GetOneUniform();

  private:
    MultipleRecursiveGenerator* mrg1;
    MultipleRecursiveGenerator* mrg2;
};


#endif
