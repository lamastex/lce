
// definitions for testing random number generators for ms


#include "config.h"
#include "test_rand.hpp"



#include <gsl/gsl_rng.h>

#include <iostream>
#include <iomanip>

using namespace std;



int main(int argc, char *argv[])
{
		randTesting1();
		randTesting2();
		randTesting3();
		
		return 0;
	
}

void randTesting1()
{
	
	int seed = 1234;
	{
		cout << "\nTry testing rand with seed " << seed << endl;
		srand(seed);
		cout << "RAND_MAX = " << RAND_MAX << endl;
		cout << rand() << endl;
		cout << rand() << endl;
		cout << rand() << endl;
		cout << rand() << endl;
		cout << rand() << endl;
	}
	{
		cout << "\nTry testing gsl_rng_rand with seed " << seed << endl;
		gsl_rng * r = gsl_rng_alloc (gsl_rng_rand);
		gsl_rng_set(r, seed);
		cout << "max = " << gsl_rng_max(r) << endl;
		cout << gsl_rng_get(r) << endl;
		cout << gsl_rng_get(r) << endl;
		cout << gsl_rng_get(r) << endl;
		cout << gsl_rng_get(r) << endl;
		cout << gsl_rng_get(r) << endl;
	}
	{
		cout << "\nTry testing gsl_rng_random_bsd with seed " << seed << endl;
		gsl_rng * r = gsl_rng_alloc (gsl_rng_random_bsd);
		gsl_rng_set(r, seed);
		cout << "max = " << gsl_rng_max(r) << endl;
		cout << gsl_rng_get(r) << endl;
		cout << gsl_rng_get(r) << endl;
		cout << gsl_rng_get(r) << endl;
		cout << gsl_rng_get(r) << endl;
		cout << gsl_rng_get(r) << endl;
	}
	{
		cout << "\nTry testing gsl_rng_random_libc5 with seed " << seed << endl;
		gsl_rng * r = gsl_rng_alloc (gsl_rng_random_libc5);
		gsl_rng_set(r, seed);
		cout << "max = " << gsl_rng_max(r) << endl;
		cout << gsl_rng_get(r) << endl;
		cout << gsl_rng_get(r) << endl;
		cout << gsl_rng_get(r) << endl;
		cout << gsl_rng_get(r) << endl;
		cout << gsl_rng_get(r) << endl;
	}
	{
		cout << "\nTry testing gsl_rng_random_glibc2 with seed " << seed << endl;
		gsl_rng * r = gsl_rng_alloc (gsl_rng_random_glibc2);
		gsl_rng_set(r, seed);
		cout << "max = " << gsl_rng_max(r) << endl;
		cout << gsl_rng_get(r) << endl;
		cout << gsl_rng_get(r) << endl;
		cout << gsl_rng_get(r) << endl;
		cout << gsl_rng_get(r) << endl;
		cout << gsl_rng_get(r) << endl;
	}
	{
		cout << "\nTry testing gsl_rng_rand48 with seed " << seed << endl;
		gsl_rng * r = gsl_rng_alloc (gsl_rng_rand48);
		gsl_rng_set(r, seed);
		cout << "max = " << gsl_rng_max(r) << endl;
		cout << gsl_rng_get(r) << endl;
		cout << gsl_rng_get(r) << endl;
		cout << gsl_rng_get(r) << endl;
		cout << gsl_rng_get(r) << endl;
		cout << gsl_rng_get(r) << endl;
	}
	
}

void randTesting2()
{
	
	int seed = 1234;
	
	{
		cout << "\nTry testing gsl_rng_random_glibc2 with seed " << seed << endl;
		cout << "and rand with seed " << seed << endl;
		srand(seed);
		cout << "RAND_MAX = " << RAND_MAX << endl;
		gsl_rng * r = gsl_rng_alloc (gsl_rng_random_glibc2);
		gsl_rng_set(r, seed);
		cout << "max = " << gsl_rng_max(r) << endl;
		cout << "gsl_rng_uniform(r) = " << gsl_rng_uniform(r) << " and rand()/(RAND_MAX+1.0) = " << rand()/(RAND_MAX + 1.0) << endl;
		cout << "gsl_rng_uniform(r) = " << gsl_rng_uniform(r) << " and rand()/(RAND_MAX+1.0) = " << rand()/(RAND_MAX + 1.0) << endl;
		cout << "gsl_rng_uniform(r) = " << gsl_rng_uniform(r) << " and rand()/(RAND_MAX+1.0) = " << rand()/(RAND_MAX + 1.0) << endl;
		cout << "gsl_rng_uniform(r) = " << gsl_rng_uniform(r) << " and rand()/(RAND_MAX+1.0) = " << rand()/(RAND_MAX + 1.0) << endl;
		cout << "gsl_rng_uniform(r) = " << gsl_rng_uniform(r) << " and rand()/(RAND_MAX+1.0) = " << rand()/(RAND_MAX + 1.0) << endl;
		cout << endl;
		cout << "gsl_rng_uniform(r) = " << gsl_rng_uniform(r) << " and rand()/(RAND_MAX+1.0) = " << rand()/(RAND_MAX + 1.0) << endl;
		cout << "gsl_rng_uniform(r) = " << gsl_rng_uniform(r) << " and rand()/(RAND_MAX+1.0) = " << rand()/(RAND_MAX + 1.0) << endl;
		cout << "gsl_rng_uniform(r) = " << gsl_rng_uniform(r) << " and rand()/(RAND_MAX+1.0) = " << rand()/(RAND_MAX + 1.0) << endl;
		cout << "gsl_rng_uniform(r) = " << gsl_rng_uniform(r) << " and rand()/(RAND_MAX+1.0) = " << rand()/(RAND_MAX + 1.0) << endl;
		cout << "gsl_rng_uniform(r) = " << gsl_rng_uniform(r) << " and rand()/(RAND_MAX+1.0) = " << rand()/(RAND_MAX + 1.0) << endl;
	}
	{
		cout << "\nTry testing gsl_rng_random_glibc2 with seed " << seed << endl;
		cout << "and rand with seed " << seed << endl;
		srand(seed);
		cout << "RAND_MAX = " << RAND_MAX << endl;
		gsl_rng * r = gsl_rng_alloc (gsl_rng_random_glibc2);
		gsl_rng_set(r, seed);
		cout << "max = " << gsl_rng_max(r) << endl;
		cout << "gsl_rng_uniform(r) = " << gsl_rng_uniform(r) << " and rand()/(RAND_MAX+1.0) = " << rand()/(RAND_MAX + 1.0) << endl;
		cout << "gsl_rng_uniform(r) = " << gsl_rng_uniform(r) << " and rand()/(RAND_MAX+1.0) = " << rand()/(RAND_MAX + 1.0) << endl;
		cout << "gsl_rng_uniform(r) = " << gsl_rng_uniform(r) << " and rand()/(RAND_MAX+1.0) = " << rand()/(RAND_MAX + 1.0) << endl;
		cout << "gsl_rng_uniform(r) = " << gsl_rng_uniform(r) << " and rand()/(RAND_MAX+1.0) = " << rand()/(RAND_MAX + 1.0) << endl;
		cout << "gsl_rng_uniform(r) = " << gsl_rng_uniform(r) << " and rand()/(RAND_MAX+1.0) = " << rand()/(RAND_MAX + 1.0) << endl;
		cout << "A couple of rand()/(RAND_MAX+1.0)" << endl;
		cout << rand()/(RAND_MAX + 1.0) << endl;
		cout << rand()/(RAND_MAX + 1.0) << endl;
		cout << "gsl_rng_uniform(r) = " << gsl_rng_uniform(r) << " and rand()/(RAND_MAX+1.0) = " << rand()/(RAND_MAX + 1.0) << endl;
		cout << "gsl_rng_uniform(r) = " << gsl_rng_uniform(r) << " and rand()/(RAND_MAX+1.0) = " << rand()/(RAND_MAX + 1.0) << endl;
		cout << "gsl_rng_uniform(r) = " << gsl_rng_uniform(r) << " and rand()/(RAND_MAX+1.0) = " << rand()/(RAND_MAX + 1.0) << endl;
		cout << "gsl_rng_uniform(r) = " << gsl_rng_uniform(r) << " and rand()/(RAND_MAX+1.0) = " << rand()/(RAND_MAX + 1.0) << endl;
		cout << "gsl_rng_uniform(r) = " << gsl_rng_uniform(r) << " and rand()/(RAND_MAX+1.0) = " << rand()/(RAND_MAX + 1.0) << endl;
	}
	{
		cout << "\nTry testing gsl_rng_random_glibc2 with seed " << seed << endl;
		cout << "and rand with seed " << seed << endl;
		srand(seed);
		cout << "RAND_MAX = " << RAND_MAX << endl;
		gsl_rng * r = gsl_rng_alloc (gsl_rng_random_glibc2);
		gsl_rng_set(r, seed);
		cout << "max = " << gsl_rng_max(r) << endl;
		cout << "gsl_rng_uniform(r) = " << gsl_rng_uniform(r) << " and rand()/(RAND_MAX+1.0) = " << rand()/(RAND_MAX + 1.0) << endl;
		cout << "gsl_rng_uniform(r) = " << gsl_rng_uniform(r) << " and rand()/(RAND_MAX+1.0) = " << rand()/(RAND_MAX + 1.0) << endl;
		cout << "gsl_rng_uniform(r) = " << gsl_rng_uniform(r) << " and rand()/(RAND_MAX+1.0) = " << rand()/(RAND_MAX + 1.0) << endl;
		cout << "gsl_rng_uniform(r) = " << gsl_rng_uniform(r) << " and rand()/(RAND_MAX+1.0) = " << rand()/(RAND_MAX + 1.0) << endl;
		cout << "gsl_rng_uniform(r) = " << gsl_rng_uniform(r) << " and rand()/(RAND_MAX+1.0) = " << rand()/(RAND_MAX + 1.0) << endl;
		cout << "A reset rand() seed to " << seed << endl;
		srand(seed);
		cout << "gsl_rng_uniform(r) = " << gsl_rng_uniform(r) << " and rand()/(RAND_MAX+1.0) = " << rand()/(RAND_MAX + 1.0) << endl;
		cout << "gsl_rng_uniform(r) = " << gsl_rng_uniform(r) << " and rand()/(RAND_MAX+1.0) = " << rand()/(RAND_MAX + 1.0) << endl;
		cout << "gsl_rng_uniform(r) = " << gsl_rng_uniform(r) << " and rand()/(RAND_MAX+1.0) = " << rand()/(RAND_MAX + 1.0) << endl;
		cout << "gsl_rng_uniform(r) = " << gsl_rng_uniform(r) << " and rand()/(RAND_MAX+1.0) = " << rand()/(RAND_MAX + 1.0) << endl;
		cout << "gsl_rng_uniform(r) = " << gsl_rng_uniform(r) << " and rand()/(RAND_MAX+1.0) = " << rand()/(RAND_MAX + 1.0) << endl;
	}
	
}

void randTesting3()
{
	
	int seed = 1234;
	
	{
		cout << "\nTry testing gsl_rng_random_glibc2 with seed " << seed << endl;
		cout << "and rand with seed " << seed << endl;
		cout << "and look for difference with simulating a positive " << seed << endl;
		srand(seed);
		cout << "RAND_MAX = " << RAND_MAX << endl;
		gsl_rng * r = gsl_rng_alloc (gsl_rng_random_glibc2);
		gsl_rng_set(r, seed);
		cout << "max = " << gsl_rng_max(r) << endl;
		
		bool okay = true;
		
		double unif_pos_glibc = 0.0;
		double unif_pos_rand = 0.0;
		
		while (okay) 
		{
			unif_pos_glibc = gsl_rng_uniform_pos(r);
			
			unif_pos_rand = rand()/(RAND_MAX+1.0);
			
			while( unif_pos_rand == 0.0 ) {
				unif_pos_rand = rand()/(RAND_MAX+1.0);
			}
					
			if (abs(unif_pos_glibc - unif_pos_rand) > 0.0000000000000001) {
				cout << setprecision(17) << endl;
				cout << "unif_pos_glibc = " << unif_pos_glibc << endl;
				cout << "unif_pos_rand = " << unif_pos_rand << endl;
				okay = false;
			}
			
		}
	}
}
