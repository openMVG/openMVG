#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <pthread.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>

#include "bitops.h"
#include "types.h"
#include "linscan.h"


/**
 * Performs kNN using linear scan in Hamming distance between codes and queries
 * Inputs:
 *   N: number of codes in the db
 *   NQ: number of queries to be answered
 *   B: number of bits in the db/query codes that should be taken into account in Hamming distance
 *   K: number of results to be returned for each query ie, k in kNN
 *   codes: an array of UINT8 storing the db codes
 *   queries: an array of UINT8 storing the query codes
 *   dim1codes: number of words in the database codes -- very likely dim1codes = B/8
 *   dim2codes: number of words in the query codes -- very likely dim2codes = B/8
 * Outputs:
 *   counter: int[(B+1)*N], stores the number of db items with different Hamming 
 *     distances (ranging from 0 to B) from each query.
 *   res: int[K*N], stores the ids of K nearest neighbors for each query
 */
void linscan_query(UINT32 *counter, UINT32 *res, UINT8 *codes, UINT8 *queries, int N, UINT32 NQ, int B, unsigned int K,
		   int dim1codes, int dim1queries) {

    int B_over_8 = B / 8;
    unsigned int *cum_counter;    // cumulative counter
    UINT32 *ind;         // stores indices arranged based on thir Hamming distances
    UINT8 *pqueries = queries;
    UINT32 *pcounter = counter;
    UINT32 *pres = res;

    memset(counter, 0, (B+1)*NQ*sizeof(*counter));

    unsigned int i=0;
#ifndef SINGLE_CORE
#pragma omp parallel shared(i) private(cum_counter, ind, pqueries, pres, pcounter)
#endif
    {
	cum_counter = new unsigned int[B+1];
	ind = new UINT32[K*(B+1)];
#ifndef SINGLE_CORE
#pragma omp for
#endif
	for (i=0; i<NQ; i++) {
#ifndef SINGLE_CORE
	    pqueries = queries + (UINT64)i*(UINT64)dim1queries;
	    pres = res + (UINT64)i*(UINT64)K;
	    pcounter = counter + (UINT64)i*(UINT64)(B+1);
#endif

	    UINT8 *pcodes = codes;
	    for (int j=0; j<N; j++, pcodes += dim1codes) {
		int h = match(pcodes, pqueries, B_over_8);
		if (h > B || h < 0) {
		    printf("Wrong Hamm distance\n");
		    exit(1);
		}
		if (pcounter[h]++ < K)
		    ind[K*h + pcounter[h] - 1] = j;
	    }
	    
	    cum_counter[0] = pcounter[0];
	    int uptoj = 0;
	    for (int j=1; j<=B; j++) {
		cum_counter[j] = cum_counter[j-1] + pcounter[j];
		if (cum_counter[j] >= K && cum_counter[j-1] < K)
		    uptoj = j;
	    }
	    
	    cum_counter[uptoj] = K;	// so we stop at K
	    
	    for (int h=0; h<=uptoj; h++) {
		int ind0 = h == 0 ? 0 : cum_counter[h-1];
		
		for (unsigned int j=ind0; j<cum_counter[h]; j++)
		    pres[j] = ind[K*h + j - ind0];
	    }
	    
#ifdef SINGLE_CORE
	    pres += K;
	    pcounter += B+1;
	    pqueries += dim1queries;
#endif
	}
	delete [] cum_counter;
	delete [] ind;
    }
}
