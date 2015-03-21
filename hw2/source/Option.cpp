
#include "Option.h"
#include "Normals.h"
#include <cmath>

double
UpAndOutCall::Price(double Spot, double Vol, double r, double d) const
{
  double pre1 = (r - d + 0.5 * Vol * Vol) * Expiry;
  double pre2 = (r - d - 0.5 * Vol * Vol) * Expiry;
  double pre0 = Vol * sqrt(Expiry);

  double d1_x_K = (log(Spot/Strike) + pre1) / pre0;
  double d2_x_K = (log(Spot/Strike) + pre2) / pre0;

  double d1_x_B = (log(Spot/Barrier) + pre1) / pre0;
  double d2_x_B = (log(Spot/Barrier) + pre2) / pre0;

  double d1_B_x = (log(Barrier/Spot) + pre1) / pre0;
  double d2_B_x = (log(Barrier/Spot) + pre2) / pre0;

  double d1_B2_Kx = (log(Barrier * Barrier / (Strike * Spot)) + pre1) / pre0;
  double d2_B2_Kx = (log(Barrier * Barrier / (Strike * Spot)) + pre2) / pre0;

  double discount = exp(-(r - d) * Expiry);
  double volPower = 2 * (r - d) / (Vol * Vol);
  return Spot * (CumulativeNormal(d1_x_K) - CumulativeNormal(d1_x_B))
     - discount * Strike * (CumulativeNormal(d2_x_K) - CumulativeNormal(d2_x_B))
     - Barrier * pow(Spot/Barrier, -volPower) * (CumulativeNormal(d1_B2_Kx) - CumulativeNormal(d1_B_x))
     + discount * Strike * pow(Spot/Barrier, 1 - volPower)
          * (CumulativeNormal(d2_B2_Kx) - CumulativeNormal(d2_B_x));
}
