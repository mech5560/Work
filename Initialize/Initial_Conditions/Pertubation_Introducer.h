#ifndef PERTUBATION_INTRODUCER_H
#define PERTUBATION_INTRODUCER_H

#include <stdlib.h>
#include <math.h>

#define pi 4.*atan(1.)



inline double drand()
{
  return (2.*((double) rand() / (RAND_MAX+1.0)) - 1.);
}


#endif
