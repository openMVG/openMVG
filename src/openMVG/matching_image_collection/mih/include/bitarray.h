#ifndef __BITARRAY_H
#define __BITARRAY_H

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "types.h"


class bitarray {

 public:

    UINT32 *arr;
    UINT32 length;
	
    bitarray()	{
	arr = NULL;
	length = 0;
    }
	
    bitarray(UINT64 _bits) {
	init(_bits);
    }
	
    void init(UINT64 _bits) {
	length = (UINT32)ceil(_bits/32.00);
	arr = new UINT32[length];
	erase();
    }
	
    ~bitarray() {
	if (arr)
	    delete[] arr;
    }
	
    inline void flip(UINT64 index) {
	arr[index >> 5] ^= ((UINT32)0x01) << (index % 32);
    }

    inline void set(UINT64 index) {
	arr[index >> 5] |= ((UINT32)0x01) << (index % 32);
    }
	
    inline UINT8 get(UINT64 index) {
	return (arr[index >> 5] & (((UINT32)0x01) << (index % 32))) != 0;
    }
	
    inline void erase() {
	memset(arr, 0, sizeof(UINT32) * length);
    }

};

#endif
