#ifndef RESULT_H__
#define RESULT_H__

#include "types.h"

struct result_t {
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
    UINT32 **res;
    UINT32 **nres;
    double **stats;
};

#endif
