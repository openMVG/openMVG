/* dynamic array of 32-bit integers
 * arr[0]   : array size (number of items inserted in the array)
 * arr[1]   : array capacity
 * arr[2..] : array content
 */

#include <stdio.h>
#include <string.h>
#include <algorithm>
#include "array32.h"

/* no need for the static keywork in the definition */
double Array32::ARRAY_RESIZE_FACTOR = 1.1;    // minimum is 1.0
double Array32::ARRAY_RESIZE_ADD_FACTOR = 1;  // minimum is 1

void Array32::set_array_resize_factor(double arf) {
  ARRAY_RESIZE_FACTOR = arf;
}

Array32::Array32() {
  arr = NULL;
}

Array32::Array32(int capacity) {
  arr = (UINT32*) calloc ((2+capacity), sizeof(UINT32));
  arr[0] = 0;
  arr[1] = capacity;
}

Array32& Array32::operator = (const Array32 &rhs) {
  if (&rhs != this)
    this->arr = rhs.arr;
}

Array32::~Array32 () {
  cleanup();
}

void Array32::cleanup () {
  free(arr);			// sometimes call free on NULL pointers, but that should be fine.
}


void Array32::expand(int newsize) {
  if (arr[1] < newsize) {
    arr[1] = newsize;
    arr[0] = newsize;
    arr = (UINT32*) realloc(arr, sizeof(UINT32)*(2 + arr[1]));
  }
}

void Array32::push(UINT32 data) {
  if (arr) {
    if (arr[0] == arr[1]) {
      arr[1] = std::max(ceil(arr[1]*ARRAY_RESIZE_FACTOR), arr[1]+ARRAY_RESIZE_ADD_FACTOR);
      arr = (UINT32*) realloc(arr, sizeof(UINT32)*(2 + arr[1]));
    }
    arr[2 + arr[0]++] = data;
  } else {
    arr = (UINT32*) malloc ((2+ARRAY_RESIZE_ADD_FACTOR)*sizeof(UINT32));
    arr[0] = 1;
    arr[1] = 1;
    arr[2] = data;
  }
}

void Array32::insert(UINT32 index, UINT32 data) {
  if (arr) {
    if (arr[0] == arr[1]) {
      arr[1] = ceil(arr[0]*ARRAY_RESIZE_FACTOR);
      arr = (UINT32*) realloc (arr, sizeof(UINT32)*(2 + arr[1]));
    }

    // if (index+1 < arr[0])
    memmove(arr+(2+index)+1, arr+(2+index), (arr[0]-index)*sizeof(*arr));
    // /* alternative to memmove */
    // for (int i=2+arr[0]; i>2+index; i--)
    //     arr[i] = arr[i-1];

    arr[2+index] = data;
    arr[0]++;

    // printf("index: %d, data: %d, array: ", index, data);
    // print();
  } else {
    arr = (UINT32*) malloc (3*sizeof(UINT32));
    arr[0] = 1;
    arr[1] = 1;
    arr[2] = data;
  }
}



UINT32* Array32::data() {
  return arr? arr + 2 : NULL;
}

UINT32 Array32::size() {
  return arr ? arr[0] : 0;
}

UINT32 Array32::capacity() {
  return arr ? arr[1] : 0;
}

void Array32::print() {
  for (int i=0; i<size(); i++) {
    printf("%d, ", arr[i+2]);
  }
  printf("\n");
}

void Array32::init(int size) {
  if (arr == NULL) {
    arr = (UINT32*) malloc ((2+size)*sizeof(UINT32));
    arr[0] = 0;
    arr[1] = size;
  }
}
