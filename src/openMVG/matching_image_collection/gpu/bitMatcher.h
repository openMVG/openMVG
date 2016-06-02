#ifndef BIT_MATCHER_H
#define BIT_MATCHER_H
void bitMatch(unsigned int*, unsigned int*, int, int, int, int*, int, cudaStream_t);
void getMatches(int, int*, int*, cudaStream_t);
#endif
