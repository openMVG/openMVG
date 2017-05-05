/*

Copyright (c) 2010--2011, Stephane Magnenat, ASL, ETHZ, Switzerland
You can contact the author at <stephane at magnenat dot net>

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the <organization> nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL ETH-ASL BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "nabo_private.h"
#include "index_heap.h"

/*!	\file brute_force_cpu.cpp
	\brief brute force search, cpu implementation
	\ingroup private
*/

namespace Nabo
{
	using namespace std;
	
	template<typename T, typename CloudType>
	BruteForceSearch<T, CloudType>::BruteForceSearch(const CloudType& cloud, const Index dim, const unsigned creationOptionFlags):
		NearestNeighbourSearch<T, CloudType>::NearestNeighbourSearch(cloud, dim, creationOptionFlags)
	{
#ifdef EIGEN3_API
		const_cast<Vector&>(this->minBound) = cloud.topRows(this->dim).rowwise().minCoeff();
		const_cast<Vector&>(this->maxBound) = cloud.topRows(this->dim).rowwise().maxCoeff();
#else // EIGEN3_API
		// compute bounds
		for (int i = 0; i < cloud.cols(); ++i)
		{
			const Vector& v(cloud.block(0,i,this->dim,1));
			const_cast<Vector&>(this->minBound) = this->minBound.cwise().min(v);
			const_cast<Vector&>(this->maxBound) = this->maxBound.cwise().max(v);
		}
#endif // EIGEN3_API
	}
	

	template<typename T, typename CloudType>
	unsigned long BruteForceSearch<T, CloudType>::knn(const Matrix& query, IndexMatrix& indices, Matrix& dists2, const Index k, const T epsilon, const unsigned optionFlags, const T maxRadius) const
	{
		const Vector maxRadii(Vector::Constant(query.cols(), maxRadius));
		return knn(query, indices, dists2, maxRadii, k, epsilon, optionFlags);
	}
	
	template<typename T, typename CloudType>
	unsigned long BruteForceSearch<T, CloudType>::knn(const Matrix& query, IndexMatrix& indices, Matrix& dists2, const Vector& maxRadii, const Index k, const T /*epsilon*/, const unsigned optionFlags) const
	{
		checkSizesKnn(query, indices, dists2, k, optionFlags, &maxRadii);
		
		const bool allowSelfMatch(optionFlags & NearestNeighbourSearch<T>::ALLOW_SELF_MATCH);
		const bool sortResults(optionFlags & NearestNeighbourSearch<T>::SORT_RESULTS);
		const bool collectStatistics(creationOptionFlags & NearestNeighbourSearch<T>::TOUCH_STATISTICS);
		
		IndexHeapSTL<Index, T> heap(k);
		
		for (int c = 0; c < query.cols(); ++c)
		{
			const T maxRadius(maxRadii[c]);
			const T maxRadius2(maxRadius * maxRadius);
			const Vector& q(query.block(0,c,dim,1));
			heap.reset();
			for (int i = 0; i < this->cloud.cols(); ++i)
			{
				const T dist(dist2<T>(this->cloud.block(0,i,dim,1), q));
				if ((dist <= maxRadius2) &&
					(dist < heap.headValue()) &&
					(allowSelfMatch || (dist > numeric_limits<T>::epsilon())))
					heap.replaceHead(i, dist);
			}
			if (sortResults)
				heap.sort();	
			heap.getData(indices.col(c), dists2.col(c));
		}
		if (collectStatistics)
			return (unsigned long)query.cols() * (unsigned long)this->cloud.cols();
		else
			return 0;
	}
	
	template struct BruteForceSearch<float>;
	template struct BruteForceSearch<double>;
	template struct BruteForceSearch<float, Eigen::Matrix3Xf>;
	template struct BruteForceSearch<double, Eigen::Matrix3Xd>;
	template struct BruteForceSearch<float, Eigen::Map<const Eigen::Matrix3Xf, Eigen::Aligned> >;
	template struct BruteForceSearch<double, Eigen::Map<const Eigen::Matrix3Xd, Eigen::Aligned> >;
}
