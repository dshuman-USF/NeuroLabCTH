//g++ -Wall -std=gnu++11 -O2 -o gen_test_data gen_test_data.cpp `pkg-config --cflags --libs gsl`
#include <stdio.h>
#include <cstdio>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <random>

const int CHANS=200;

double probs[CHANS];

int
main (int argc, char **argv)
{
  gsl_rng_env_setup ();
  gsl_rng *rng = gsl_rng_alloc (gsl_rng_default);
   
  std::random_device rd; // obtain a random number from hardware
  std::mt19937 eng(rd()); // seed the generator
  std::uniform_int_distribution<> distr(350000, 400000); // define the range

  probs[0] = 425199. / 111588909; // original
  probs[1] = 374403. / 111588909;

  for (int prb=2; prb<CHANS*2; prb += 2)
  {
     probs[prb] = distr(eng) / 111588909.0;
     probs[prb+1] = distr(eng) / 111588909.0;
  }

  int chan;
  int max_tick = 10000000;
  int cycle_ticks = 1075 * 50;
  int e_phase = 1075 * 36;
//  double p1 = 425199. / 111588909;
//  double p2 = 374403. / 111588909;

  printf ("%5d%10d\n", 33, 3333333);
  printf ("%5d%10d\n", 33, 3333333);
  printf ("%5d%10d\n", 58, 0);  // start of period
  int t;
  for (t = 0; t < max_tick; t++) 
  {
    if (t % cycle_ticks == e_phase)
      printf ("%5d%10d\n", 97, t);
    if (t % cycle_ticks == 0)
      printf ("%5d%10d\n", 98, t);
    for (chan=101 ; chan < CHANS; chan+= 2)
    { 
        // errr. . . do I need a lot of p values?
//      if (gsl_rng_uniform (rng) < p1)
//        printf ("%5d%10d\n", chan, t);
//      if (gsl_rng_uniform (rng) < p2)
//        printf ("%5d%10d\n", chan+1, t);
      if (gsl_rng_uniform (rng) < probs[chan])
        printf ("%5d%10d\n", chan, t);
      if (gsl_rng_uniform (rng) < probs[chan+1])
        printf ("%5d%10d\n", chan+1, t);
    }
  }
  printf ("%5d%10d\n", 58, t);  // end of period

  return 0;
}

