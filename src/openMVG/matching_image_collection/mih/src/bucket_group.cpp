#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "types.h"
#include "bitops.h"
#include "array32.h"
#include "bucket_group.h"

BucketGroup::BucketGroup() {
    empty = 0;
    group = NULL;
}

BucketGroup::~BucketGroup() {
    if (group != NULL)
	delete group;
}

void BucketGroup::count_insert(int subindex, UINT32 data) {
    // printf("%ld\n", (long int)group);
    // printf("(%d, %d)\n", subindex, data);

    if (group == NULL) {
	group = new Array32();
	group->push(11111);
    }

    // group->push(subindex);

    UINT32 lowerbits = ((UINT32)1 << subindex) - 1;
    int howmany = popcnt(empty & lowerbits);
    if (!(empty & ((UINT32)1 << subindex))) {
    	group->insert(howmany+1, 1);
    	empty |= (UINT32)1 << subindex;
    } else
    	group->arr[2 + howmany+1]++;
    // printf("howmany: %d, empty: %d\n", howmany, empty);
}

void BucketGroup::data_insert(int subindex, UINT32 data) {
    // sanity check
    if (group == NULL) {
	printf("something is prob wrong, group is NULL.\n");
	exit(1);
    }

    if (group->arr[2] == 11111) {
	group->arr[2] = 0;
	allocate_mem_based_on_counts();
    }

    // group is certainly not NULL
    int totones = popcnt(empty);
    UINT32 lowerbits = ((UINT32)1 << subindex) - 1;
    int howmany = popcnt(empty & lowerbits);

    group->arr[2 + totones+1 + group->arr[2+howmany+1]++] = data;
}

void BucketGroup::allocate_mem_based_on_counts() {
    if (group == NULL)
	return;

    // int n = group->size()-1; // -1 because of the number pushed at the beginning (11111 or 0)
    // int count_subindex[32+1];
    // memset(count_subindex, 0, (32+1)*sizeof(int));

    // for (int i=1; i<group->size(); i++) {
    // 	count_subindex[group->arr[2+i]+1]++;
    // 	empty |= (UINT32)1 << group->arr[2+i];
    // }
    // for (int i=1; i<33; i++)
    // 	count_subindex[i] += count_subindex[i-1];

    // int c = 0;
    // for (int i=0; i<32; i++)
    // 	if (count_subindex[i] != count_subindex[i+1])
    // 	    group->arr[2 + 1+c++] = count_subindex[i];
    // group->expand(n + c+1);

    // //

    int totones = popcnt(empty);
    assert(totones+1 == group->size());

    for (int i=0; i<totones; i++)
    	group->arr[2 + i+1] += group->arr[2 + i];

    group->expand(group->arr[2 + totones] + totones+1);

    for (int i=totones-1; i>=0; i--)
    	group->arr[2 + i+1] = group->arr[2 + i];    
}


void BucketGroup::lazy_insert(int subindex, UINT32 data) {
    if (group == NULL) {
	group = new Array32();
    }

    empty |= (UINT32)1 << subindex;
    group->push(data);
}

void BucketGroup::cleanup_insert(UINT8* dataset, int m, int k, int mplus, int b, int dim1codes) {
    int totones = popcnt(empty);
    Array32 *group2 = new Array32(group->size() + totones+1);
    int n = group->size();
    int *subindices = new int[n];
    UINT64 *chunks = new UINT64[m];
    int cumcount[33];
    memset(cumcount, 0, 33*sizeof(int));

    for (int i=0; i<n; i++) {
	UINT8 *current = dataset+(size_t)(group->arr[2+i])*(size_t)dim1codes;
	split(chunks, current, m, mplus, b);
	subindices[i] = (int)(chunks[k] & 31);
	//	printf("%d: %d\n", i, subindices[i]);
	cumcount[1 + subindices[i]]++;
    }

    for (int i=0; i<32; i++)
	cumcount[i+1] += cumcount[i];

    // int howmany[32];
    // int lowerbits = 0;
    // for (int i=0; i<32; i++) {
    // 	howmany[i] = popcnt(empty & lowerbits);
    // 	lowerbits |= 1 << i;
    // }

    group2->arr[0] = group2->arr[1];
    for (int i=0; i<n; i++)
	group2->arr[2 + totones+1 + cumcount[subindices[i]]++] = group->arr[2+i];

    // group2->print();

    group2->arr[2] = 0;
    int c = 0;
    for (int i=0; i<32; i++)
	if (empty & (1<<i))
	    group2->arr[2+1+c++] = cumcount[i];
    
    free(group);
    group = group2;
    delete[] chunks;
    delete[] subindices;
}

void BucketGroup::insert(int subindex, UINT32 data) {
    if (group == NULL) {
	group = new Array32();
	group->push(0);
    }

    /* arr[2 .. totones] stores the indices for the position of the corresponding sub-buckets in arr */

    UINT32 lowerbits = ((UINT32)1 << subindex) - 1;
    int howmany = popcnt(empty & lowerbits);
    if (!(empty & ((UINT32)1 << subindex))) {
	group->insert(howmany, group->arr[howmany+2]);
	empty |= (UINT32)1 << subindex;
    }

    /* the number of ones in empty (after a potential insert above):
       the number of sub-buckets that are non-empty */
    int totones = popcnt(empty);

    group->insert(totones+1+group->arr[2+howmany+1], data);
    for (int i=howmany+1; i<totones+1; i++)
	group->arr[2+i]++;
}

UINT32* BucketGroup::query(int subindex, int *size) {
    if (empty & ((UINT32)1 << subindex)) {
	UINT32 lowerbits = ((UINT32)1 << subindex) - 1;
	int howmany = popcnt(empty & lowerbits);
	int totones = popcnt(empty);
	*size = group->arr[2+howmany+1]-group->arr[2+howmany];
	return group->arr + 2 + totones+1 + group->arr[2+howmany];
    } else {
	*size = 0;
	return NULL;
    }
}
