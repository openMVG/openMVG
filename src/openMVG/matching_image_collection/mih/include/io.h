#ifndef IO_H__
#define IO_H__

#include "result.h"

const int MAXQ = 10000;

/* 
 *  Uses hdf5 library to load a dataset (varStr) from a file (filename).
 */
void load_bin_codes(const char *filename, const char *varStr, UINT8 *var, 
                    UINT32 *N, int *B, int start0=0);

void saveVar(const char *filename, const char *varStr, UINT8 *var, int N, int B,
             int overwrite = 1);

void saveRes(const char *filename, const char *varStr, const result_t *result,
             int n, int overwrite = 1);

#endif
