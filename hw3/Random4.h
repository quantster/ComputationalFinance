// Nitish Kanabar

#ifndef RANDOM4_H
#define RANDOM4_H

#include "Arrays.h"

// UNIFORM_GENERATOR CODES:
// 1 -> ParkMiller Uniforms
// 2 -> LEcuyer Uniforms
// else -> RandomRand Uniforms
#define UNIFORM_GENERATOR 0

#ifndef UNIFORM_GENERATOR
#define UNIFORM_GENERATOR 0
#endif

#if UNIFORM_GENERATOR == 1
#include "ParkMiller.h"
#elif UNIFORM_GENERATOR == 2
#include "LEcuyer.h"
#else
#include "Random3.h"
#endif


#if UNIFORM_GENERATOR == 1
class FishmanGenerator : public RandomParkMiller
#elif UNIFORM_GENERATOR == 2
class FishmanGenerator : public RandomLEcuyer
#else
class FishmanGenerator : public RandomRand
#endif
{
  public:
    FishmanGenerator(unsigned long Dimensionality);
    virtual void GetGaussians(MJArray& variates);
    virtual double GetOneGaussian();
};

#if UNIFORM_GENERATOR == 1
class BoxMullerGenerator : public RandomParkMiller
#elif UNIFORM_GENERATOR == 2
class BoxMullerGenerator : public RandomLEcuyer
#else
class BoxMullerGenerator : public RandomRand
#endif
{
  public:
    BoxMullerGenerator(unsigned long Dimensionality);
    virtual void GetGaussians(MJArray& variates);
    virtual double GetOneGaussian();
};


#if UNIFORM_GENERATOR == 1
class InverseTransformGenerator : public RandomParkMiller
#elif UNIFORM_GENERATOR == 2
class InverseTransformGenerator : public RandomLEcuyer
#else
class InverseTransformGenerator : public RandomRand
#endif
{
  public:
    InverseTransformGenerator(unsigned long Dimensionality);
    virtual void GetGaussians(MJArray& variates);
    virtual double GetOneGaussian();
};

#endif
