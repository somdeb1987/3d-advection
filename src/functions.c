#include <math.h>

double absolute(double x)
{
  if (x < 0) {
    return -x;
  } else {
    return x;
  }
}

double raiseto(double x, double a)
{
  return exp(a*log(x));
}
