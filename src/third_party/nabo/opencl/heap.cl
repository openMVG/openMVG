T heapHeadValue(const HeapEntry* heap)
{
	return heap->value;
}

void heapHeadReplace(HeapEntry* heap, const int index, const T value, const uint K)
{
	uint i = 0;
	for (; i < K - 1; ++i)
	{
		if (heap[i + 1].value > value)
			heap[i] = heap[i + 1];
		else
			break;
	}
	heap[i].value = value;
	heap[i].index = index;
}

void heapInit(HeapEntry* heap, const uint K)
{
	for (uint i = 0; i < K; ++i)
		heap[i].value = HUGE_VALF;
}

void heapCopy(global int* indices, global T* dists2, const HeapEntry* heap, const uint K)
{
	for (uint i = 0; i < K; ++i)
	{
		*indices++ = heap[K-i-1].index;
		*dists2++ = heap[K-i-1].value;
	}
}

void heapSort(const HeapEntry* heap)
{
}
