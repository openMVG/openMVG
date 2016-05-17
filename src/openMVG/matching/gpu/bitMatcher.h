#ifndef BIT_MATCHER_H
#define BIT_MATCHER_H
void bitMatcher(unsigned int*, unsigned int*, int, int, int, int*, int, cudaStream_t, cudaEvent_t);
void getMatches(int, int*, int*);
#endif
