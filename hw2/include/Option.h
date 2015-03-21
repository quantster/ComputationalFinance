/*
	BarrierOption.h
*/

#ifndef BARRIER_OPTION_H_
#define BARRIER_OPTION_H_

#include "minmax.h"

class VanillaCall
{
  public:
    VanillaCall(double Strike_, double Expiry_)
      : Strike(Strike_), Expiry(Expiry_) {}

    double GetExpiry() const
    {
      return Expiry;
    }

    void Reset() {}
    bool IsKnockedOut(double) const { return false; }
    void Record(double) {}

    double GetStrike() const
    {
      return Strike;
    }

    double Payoff(double Spot) const
    {
      return max(Spot - Strike, 0.0);
    }

  private:
    double Strike;
    double Expiry;
};


class UpAndOutCall
{
	public:
		UpAndOutCall(double Strike_, double Barrier_, double Expiry_)
			: Strike(Strike_), Barrier(Barrier_), Expiry(Expiry_), isOut(false) {}
	 
	  void Record(double Spot)
		{
			isOut = isOut || (Spot > Barrier);
		}

    double GetStrike() const
    {
      return Strike;
    }

		double Payoff(double Spot) const
		{
			return isOut ? 0.0 : max(Spot - Strike, 0.0);
		}

    double GetExpiry() const
    {
      return Expiry;
    }

    bool IsKnockedOut(double Spot)
    {
			isOut = isOut || (Spot > Barrier);
      return isOut;
    }

    void Reset()
    {
      isOut = false;
    }


    double  Price(double Spot, double Vol, double r, double d) const;

	private:
		double Strike;
		double Barrier;
    double Expiry;
		bool isOut;
};

#endif
