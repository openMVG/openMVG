/* -*- mode: C++; indent-tabs-mode: nil; -*-
 *
 * This file is a part of LEMON, a generic C++ optimization library.
 *
 * Copyright (C) 2003-2009
 * Egervary Jeno Kombinatorikus Optimalizalasi Kutatocsoport
 * (Egervary Research Group on Combinatorial Optimization, EGRES).
 *
 * Permission to use, modify and distribute this software is granted
 * provided that this copyright notice appears in all copies. For
 * precise terms see the accompanying LICENSE file.
 *
 * This software is provided "AS IS" with no warranty of any kind,
 * express or implied, and with no claim as to its suitability for any
 * purpose.
 *
 */

#include <lemon/random.h>
#include "test_tools.h"

int seed_array[] = {1, 2};

int rnd_seq32[] = {
2732, 43567, 42613, 52416, 45891, 21243, 30403, 32103, 
62501, 33003, 12172, 5192, 32511, 50057, 43723, 7813, 
23720, 35343, 6637, 30280, 44566, 31019, 18898, 33867, 
5994, 1688, 11513, 59011, 48056, 25544, 39168, 25365, 
17530, 8366, 27063, 49861, 55169, 63848, 11863, 49608
};
int rnd_seq64[] = {
56382, 63883, 59577, 64750, 9644, 59886, 57647, 18152, 
28520, 64078, 17818, 49294, 26424, 26697, 53684, 19209, 
35404, 12121, 12837, 11827, 32156, 58333, 62553, 7907, 
64427, 39399, 21971, 48789, 46981, 15716, 53335, 65256, 
12999, 15308, 10906, 42162, 47587, 43006, 53921, 18716
};

void seq_test() {
  for(int i=0;i<5;i++) {
    lemon::Random32 r32(i);
    lemon::Random64 r64(i);
    for(int j=0;j<8;j++) {
      check(r32[65536]==rnd_seq32[i*8+j], "Wrong random sequence");
      check(r64[65536]==rnd_seq64[i*8+j], "Wrong random sequence");
    }
  }
}


int main()
{
  double a=lemon::rnd();
  check(a<1.0&&a>0.0,"This should be in [0,1)");
  a=lemon::rnd.gauss();
  a=lemon::rnd.gamma(3.45,0);
  a=lemon::rnd.gamma(4);
  //Does gamma work with integer k?
  a=lemon::rnd.gamma(4.0,0);
  a=lemon::rnd.poisson(.5);

  lemon::rnd.seed(100);
  lemon::rnd.seed(seed_array, seed_array +
                  (sizeof(seed_array) / sizeof(seed_array[0])));

  seq_test();
  return 0;
}
