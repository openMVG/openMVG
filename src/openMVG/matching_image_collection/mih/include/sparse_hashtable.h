#ifndef SPHASHTABLE_H__
#define SPHASHTABLE_H__

#include <stdio.h>
#include <math.h>
#include "types.h"
#include "bucket_group.h"

class SparseHashtable {

 private:
    static const int MAX_B;	// Maximum bits per key before folding the table	

    BucketGroup *table;		// Bins (each bin is an Array object for duplicates of the same key)

 public:

    int b;			// Bits per index
	
    UINT64 size;		// Number of bins

    SparseHashtable();

    ~SparseHashtable();

    int init(int _b);
	
    void insert(UINT64 index, UINT32 data);

    void lazy_insert(UINT64 index, UINT32 data);

    UINT32* query(UINT64 index, int* size);

    void cleanup_insert(UINT8* dataset, int m, int k, int mplus, int b, int dim1codes);

    void count_insert(UINT64 index, UINT32 data);

    void allocate_mem_based_on_counts();

    void data_insert(UINT64 index, UINT32 data);

};

#endif
