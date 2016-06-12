#ifndef RESULTH5_H__
#define RESULTH5_H__

#include "hdf5.h"
#include "types.h"

/* The sequence of the fields of resulth5_t and their types should be
   exactly the same as result_t. The pointer types change to hdf5 pointer
   objects. */
struct resulth5_t {
    int n;
    int nq;
    int k;
    int b;
    int m;
    int q0;
    int q1;
    double wt;
    double cput;
    double vm;
    double rss;
    hobj_ref_t res;
    hobj_ref_t nres;
    hobj_ref_t stats;
};

#endif
