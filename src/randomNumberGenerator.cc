//&&&&&&&&&&&&&&&&
//& Header Files &
//&&&&&&&&&&&&&&&&
#include "randomNumberGenerator.hh"

//&&&&&&&&&&&&&&&&&&&&&&&&
//& Function Definitions &
//&&&&&&&&&&&&&&&&&&&&&&&&

void rngSeed(int generatorFlag, int seed, gsl_rng * r)
{
//This function seeds the rng()

	long unsigned int seedUL;

	if (generatorFlag==0)
	{	
		srand(seed);
	}
	else
	{
		seedUL = (long unsigned int) seed;
		gsl_rng_set(r, seedUL);
	}

	return;
}

double rng(int generatorFlag, gsl_rng * r)
{
//This function will return some random numbers

	if (generatorFlag==0)		//cstdlib rand()
	{	
		return  ((double) rand())/( (double) RAND_MAX + 1.0);
	}
	else
	{	
		return gsl_rng_uniform(r);
	}
}
