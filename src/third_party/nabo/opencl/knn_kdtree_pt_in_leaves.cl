#define INVALID_NODE -1

#define OP_BEGIN_FUNCTION 0
#define OP_REC1 1
#define OP_REC2 2

typedef struct { 
	int dim; // -1 == invalid, <= -2 = index of pt
	T cutVal; //!< for split node, split value 
} Node;

typedef struct {
	uint op;
	size_t n;
	size_t other_n;
	T rd;
	T old_off;
	T new_off;
} StackEntry;

size_t childLeft(const size_t pos) { return 2*pos + 1; }
size_t childRight(const size_t pos) { return 2*pos + 2; }

void offInit(T *off)
{
	for (uint i = 0; i < DIM_COUNT; ++i)
		off[i] = 0;
}

// for cloud and result, use DirectAccessBit and stride
// preconditions:
// 		K < MAX_K
//		stack_ptr < MAX_STACK_DEPTH
kernel void knnKDTree(	const global T* cloud,
						const global T* query,
						global int* indices,
						global T* dists2,
						const uint K,
						const T maxError,
						const T maxRadius2,
						const uint optionFlags,
						const uint indexStride,
						const uint dists2Stride,
						const uint pointCount,
#ifdef TOUCH_STATISTICS
						global uint* touchStatistics,
#endif
						const global Node* nodes
 					)
{
	StackEntry stack[MAX_STACK_DEPTH];
	HeapEntry heap[MAX_K];
	T off[DIM_COUNT];

#ifdef TOUCH_STATISTICS
	uint visitCount = 0;
#endif

	const size_t queryId = get_global_id(0);
	const bool allowSelfMatch = optionFlags & ALLOW_SELF_MATCH;
	const bool doSort = optionFlags & SORT_RESULTS;
	const global T* q = &query[queryId * POINT_STRIDE];
	
	heapInit(heap, K);
	offInit(off);
	
	uint stackPtr = 1;
	stack[0].op = OP_BEGIN_FUNCTION;
	stack[0].n = 0;
	stack[0].rd = 0;

	while (stackPtr != 0)
	{
		--stackPtr;
		StackEntry* s = stack + stackPtr;
		const size_t n = s->n;
		const global Node* node = nodes + n;
		const int cd = node->dim;
		switch (stack[stackPtr].op)
		{
			case OP_BEGIN_FUNCTION:
			if (cd < 0)
			{
				if (cd != -1)
				{
					const int index = -(cd + 2);
					const global T* p = &cloud[index * POINT_STRIDE];
					T dist = 0;
					for (uint i = 0; i < DIM_COUNT; ++i)
					{
						const T diff = q[i] - p[i];
						dist += diff * diff;
					}
					if ((dist <= maxRadius2) &&
						(dist < heapHeadValue(heap) &&
						(allowSelfMatch || (dist > (T)EPSILON))))
						heapHeadReplace(heap, index, dist, K);
#ifdef TOUCH_STATISTICS
					++visitCount;
#endif
				}
			}
			else
			{
				s->old_off = off[cd];
				s->new_off = q[cd] - node->cutVal;
				(s+1)->op = OP_BEGIN_FUNCTION;
				s->op = OP_REC1;
				if (s->new_off > 0)
				{
					(s+1)->n = childRight(n);
					s->other_n = childLeft(n);
				}
				else
				{
					(s+1)->n = childLeft(n);
					s->other_n = childRight(n);
				}
				(s+1)->rd = s->rd;
				stackPtr+=2;
			}
			break;
			
			case OP_REC1:
			{
				T rdE;
				s->rd += - (s->old_off*s->old_off) + (s->new_off*s->new_off);
				if ((s->rd <= maxRadius2) &&
					(s->rd * maxError < heapHeadValue(heap)))
				{
					off[cd] = s->new_off;
					(s+1)->op = OP_BEGIN_FUNCTION;
					(s+1)->n = s->other_n;
					(s+1)->rd = s->rd;
					s->op = OP_REC2;
					stackPtr+=2;
				}
			}
			break;
			
			case OP_REC2:
			off[cd] = s->old_off;
			break;
		}
	}
	
	if (doSort)
		heapSort(heap);
	heapCopy(&indices[queryId * indexStride], &dists2[queryId * dists2Stride], heap, K);
#ifdef TOUCH_STATISTICS
	touchStatistics[queryId] = visitCount;
#endif
}

