#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <pthread.h>
#include <time.h>
#include <string.h>
#include <vector>

#include "bitops.h"
#include "types.h"
#include "linscan.h"

using namespace std;

int minmax (int *C, vector<int> cbits, int * freebits, int d) {
	int bestbit = -1;
	int min = C[0]+1; // this is the max correlation possible
	
	for (int j=0; j<d; j++) {
		if (freebits[j]) {
			int max = 0;
			for (int i=0; i<cbits.size(); i++) {
				int v = abs(C[j*d + cbits[i]]);
				if (v > max) max = v;
			}
			if (max < min) {
				min = max;
				bestbit = j;
			}
		}
	}
	
	return bestbit;	
}

int maxmax (int *C, vector<int> cbits, int * freebits, int d) {
	int bestbit = -1;
	int maxmax = 0;
	
	for (int j=0; j<d; j++) {
		if (freebits[j]) {
			int max = 0;
			for (int i=0; i<cbits.size(); i++) {
				int v = abs(C[j*d + cbits[i]]);
				if (v > max) max = v;
			}
			if (max > maxmax) {
				maxmax = max;
				bestbit = j;
			}
		}
	}
	
	return bestbit;
}

void greedyorder (int * order, UINT8 * B, size_t N, int d, int m) {
// B is database points, column major
// N is the number of points to use
// d is the number of bits to use
// m is the number of chunks to use
	
	int *C = new int [d * d];	// Correlation matrix, column major

	UINT8 mask [] = {0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80};
	UINT8 * BB = new UINT8 [N * d]; // Unpacked bits, column major
	

	for (size_t i=0; i<N; i++) {
		for (int b1=0; b1<d; b1+=8) {
			for (int b2=0; b2<8; b2++) {
				BB[i*d + b1 + b2] = (B[i*d/8 + b1/8] & mask[b2])>(UINT8)0;
			}
		}
	}
	for (int j=0; j<d; j++) 
		for (int k=0; k<d; k++) 
			C[j*d + k] = 0;
			
	int j;
	#pragma omp parallel shared(j)
	{
	#pragma omp for
		for (j=0; j<d; j++) {
			for (size_t i=0; i<N; i++) {
				for (int k=j+1; k<d; k++) {
					C[j*d + k] -= BB[i*d + j] ^ BB[i*d + k];
				}
			}
		}
	}	
	
	for (int j=0; j<d; j++)
		for (int k=0; k<=j; k++) {
			C[k*d + j] *= 2;
			C[k*d + j] += N;
			C[j*d + k] = C[k*d + j];
			if (j==k) C[j*d + k] = N;
		}
/*	
	for (int i=0; i<N; i++) {
		for (int by=0; by<d/8; by++) printf ("%5d", B[i*d/8 + by]);
		printf ("\n");
	}
	
	printf ("\n-------------------------------\n");
	
	for (int i=0; i<N; i++) {
		for (int b=0; b<d; b++) {
			if (b%8==0) printf (" ");
			printf ("%1d", BB[i*d + b]);
		}
		printf ("\n");
	}
	
	
	printf ("\n--------------------------------\n");
	
	for (int j=0; j<16; j++) {
		for (int k=0; k<16; k++) {
			printf ("%6.2f", (float)C[j*d + k]/N);
		}
		printf("\n");
	}
*/

	int* freebits = new int [d];
	for (int i=0; i<d; i++) freebits[i] = 1;
	
	vector<int>* chunks = new vector<int> [m];
	
	int bestbit;
	for (int b=0; b<d; b++) {
		int c = b % m;
		if (chunks[c].size() == 0) {
			if (c==0) bestbit = 0;
			else bestbit = maxmax(C, chunks[c-1], freebits, d);
		} else {
			bestbit = minmax(C, chunks[c], freebits, d);
		}
		chunks[c].push_back(bestbit);
		freebits[bestbit] = 0;
	}
	
	int b=0;
	for (int c=0; c<m; c++) {
		for (int i=0; i<chunks[c].size(); i++)
			order[b++] = chunks[c][i];
	}	
	
	for (int c=0; c<m; c++) {
		for (int i=0; i<chunks[c].size(); i++)
			printf("%4d", chunks[c][i]);
		printf("\n");
	}
	
	delete [] chunks;
	delete freebits;
	delete C;
	delete BB;
	
}

void reorder(UINT8 * out, UINT8 * in, size_t N, int d, int * order) {
// reorder bits according to order
// N codes, d bits each
// out and in are the same size, both are d x N, column major

	UINT8 mask [] = {0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80};
	
	memset(out, 0, (size_t)N * d/8);
	
	int * order_byte = new int [d];
	int * order_bit = new int [d];
	for (int i=0; i<d; i++) {
		order_byte[i] = order[i]/8;
		order_bit[i] = order[i]%8;
	}
	size_t i;
	#pragma omp parallel shared(i)
	{
	#pragma omp for
		for (i =0; i<N; i++) {
			for (int b=0; b<d; b++) {
				int by = b/8;
				int bi = b%8; 
				if (in[i*d/8 + order_byte[b]] & mask[order_bit[b]]) // Bit-shifting is about 10% faster, instead of mask lookup
					out[i*d/8 + by] |= mask[bi];
			}
		}
	}	
/*
	for (int i=0; i<16; i++) printf ("%4d",order[i]); printf("\n");
	for (int i=0; i<16; i++) printf ("%4d",order_byte[i]); printf("\n");
	for (int i=0; i<16; i++) printf ("%4d",order_bit[i]); printf("\n");
	printf("\n------------------------------\n");
	for (size_t i=N-5; i<N; i++) {
		for (int by=0; by<d/8; by++) {
			printf("%4d", in[i*d/8 + by]);
		}
		printf("\n");
	}
	printf("\n------------------------------\n");
	for (size_t i=N-5; i<N; i++) {
		for (int b=0; b<d; b++) {
			if (b%8 == 0) printf(" ");
			printf("%1d", (in[i*d/8 + b/8] & mask[b%8])>0 );
		}
		printf("\n");
	}
	printf("\n------------------------------\n");
	for (size_t i=N-5; i<N; i++) {
		for (int b=0; b<d; b++) {
			if (b%8 == 0) printf(" ");
			printf("%1d", (out[i*d/8 + b/8] & mask[b%8])>0 );
		}
		printf("\n");
	}
*/
	delete [] order_byte;
	delete [] order_bit;
}
