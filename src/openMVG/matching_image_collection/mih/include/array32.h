/* dynamic array of 32-bit integers
 * arr[0]   : array size
 * arr[1]   : array capacity
 * arr[2..] : array content */

#ifndef ARRAY32_H__
#define ARRAY32_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "types.h"

class Array32 {

 private:
    static double ARRAY_RESIZE_FACTOR;
    static double ARRAY_RESIZE_ADD_FACTOR;

 public:

    static void set_array_resize_factor(double arf);

    UINT32 *arr;

    Array32();

    /* initializes to zero too */
    Array32(int capacity);

    ~Array32();

    void cleanup();

    void push(UINT32 data);

    void insert(UINT32 index, UINT32 data);

    UINT32* data();

    UINT32 size();

    UINT32 capacity();

    Array32& operator= (const Array32&);	

    void print();

    void init(int size);

    void expand(int newsize);

};

#endif
