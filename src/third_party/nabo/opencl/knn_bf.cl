kernel void knnBruteForce(const global T* cloud,
						const global T* query,
						global int* indices,
						global T* dists2,
						const uint K,
						const T maxError,
						const T maxRadius2,
						const uint optionFlags,
						const uint indexStride,
						const uint dists2Stride,
						const uint pointCount
#ifdef TOUCH_STATISTICS
						,
						global uint* touchStatistics
#endif
						)
{
	HeapEntry heap[MAX_K];
	heapInit(heap, K);
	
	const size_t queryId = get_global_id(0);
	const bool allowSelfMatch = optionFlags & ALLOW_SELF_MATCH;
	const bool doSort = optionFlags & SORT_RESULTS;
	const global T* q = &query[queryId * POINT_STRIDE];
	
	for (uint index = 0; index < pointCount; ++index)
	{
		const global T* p = &cloud[index * POINT_STRIDE];
		T dist = 0;
		for (uint i = 0; i < DIM_COUNT; ++i)
		{
			const T diff = q[i] - p[i];
			dist += diff * diff;
		}
		if ((dist <= maxRadius2) &&
			(dist < heapHeadValue(heap) &&
			(allowSelfMatch || (dist > 0)))
			heapHeadReplace(heap, index, dist, K);
	}
	
	if (doSort)
		heapSort(heap);
	heapCopy(&indices[queryId * indexStride], &dists2[queryId * dists2Stride], heap, K);
	#ifdef TOUCH_STATISTICS
	touchStatistics[queryId] = pointCount;
	#endif
}
