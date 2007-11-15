#include "math.h"
#include "random.h"

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836

/* constructor & destructor */

Random::Random()
{
  seed = 123456559;
}

/* Park/Miller generator */
/*set seed*/
void Random::set_seed(int inseed)
{
  seed = inseed;
}
double Random::uniform()
{
  int k;
  double ans;

  k = seed/IQ;
  seed = IA*(seed-k*IQ) - IR*k;
  if (seed < 0) seed += IM;
  ans = AM*seed;
  return ans;
}

/* gaussian ramdom number generator */

double Random::gaussian(double sigma, double mu)
{
  int k;
  double ans;
  double sign;

  ans =  sqrt(-sigma*sigma*log(uniform()));
  sign = cos(2.0*3.141592*uniform());
  ans = mu + sign * ans;

  return ans;
}

/* exponential random number generator */

double Random::exponential(double pref)
{
  double ans;
  double sign;

  ans = -log(uniform())/pref;

  return ans;
}

