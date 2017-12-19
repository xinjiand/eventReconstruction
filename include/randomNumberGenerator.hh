#ifndef RANDOM_NUMBER_GENERATOR_H
#define RANDOM_NUMBER_GENERATOR_H

//&&&&&&&&&&&&&&&&&
//& C headerfiles &
//&&&&&&&&&&&&&&&&&
#include <stdlib.h>

//&&&&&&&&&&&&&&&&&&&&&&&&&&
//& GNU Scientific Library &
//&&&&&&&&&&&&&&&&&&&&&&&&&&
#include <gsl/gsl_rng.h>

//&&&&&&&&&&&&&&&&&&&&&&&&&
//& Function Declarations &
//&&&&&&&&&&&&&&&&&&&&&&&&&

void		rngSeed(int generatorFlag, int seed, gsl_rng * r);
double	rng(int generatorFlag, gsl_rng * r);

#endif
