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

#ifndef __INDEX_HEAP_H
#define __INDEX_HEAP_H

#include "nabo.h"
#include <vector>
#include <algorithm>
#include <limits>
#include <iostream>

/*!	\file index_heap.h
	\brief implementation of index heaps
	\ingroup private
*/

namespace Nabo
{
	//! balanced-tree implementation of heap
	/** It uses a binary heap, which provides replacement in O(log(n)),
	 * 	however the constant overhead is significative. */
	template<typename IT, typename VT>
	struct IndexHeapSTL
	{
		//! type of an index
		typedef IT Index;
		//! type of a value
		typedef VT Value;
		
		//! an entry of the heap tree
		struct Entry
		{
			IT index; //!< index of point
			VT value; //!< distance for this point
			
			//! create a new entry
			Entry(const IT index, const VT value): index(index), value(value) {}
			//! return true if e0 is of lower value than e1, false otherwise
			friend bool operator<(const Entry& e0, const Entry& e1) { return e0.value < e1.value; }
		};
		//! vector of entry, type for the storage of the tree
		typedef std::vector<Entry> Entries;
		//! vector of indices
		typedef typename Eigen::Matrix<Index, Eigen::Dynamic, 1> IndexVector;
		//! vector of values
		typedef typename Eigen::Matrix<Value, Eigen::Dynamic, 1> ValueVector;
		
		//! storage for the tree
		Entries data;
		//! number of neighbours requested
		const size_t nbNeighbours;

		
		//! Constructor
		/*! \param size number of elements in the heap */
		IndexHeapSTL(const size_t size):
			data(1, Entry(0, std::numeric_limits<VT>::infinity())),
			nbNeighbours(size)
		{
			data.reserve(size);
		}
		
		//! reset to the empty heap
		inline void reset()
		{
			data.clear();
			data.push_back(Entry(0, std::numeric_limits<VT>::infinity()));
		}
		
		//! get the largest value of the heap
		/** \return the largest value in the heap */
		inline const VT& headValue() const { return data.front().value; }
		
		//! put value into heap, replace the largest value if full
		/** \param index new point index
		 * 	\param value new distance value */
		inline void replaceHead(const Index index, const Value value)
		{

			if (data.size() == nbNeighbours)
			{	// we have enough neighbours to discard largest
				pop_heap(data.begin(), data.end());
				data.back() = Entry(index, value);
			}
			else
			{	// missing neighbours
				data.push_back(Entry(index, value));
			}
			// ensure heap
			push_heap(data.begin(), data.end());
		}
		
		//! sort the entries, from the smallest to the largest
		inline void sort()
		{
			sort_heap (data.begin(), data.end());
		}
		
		//! get the data from the heap
		/** \param indices index vector
		 * 	\param values value vector */
		template<typename DI, typename DV>
		inline void getData(const Eigen::MatrixBase<DI>& indices, const Eigen::MatrixBase<DV> & values) const
		{
			// note: we must implement this hack because of problem with reference to temporary
			// C++0x will solve this with rvalue
			// see: http://eigen.tuxfamily.org/dox-devel/TopicFunctionTakingEigenTypes.html
			// for more informations
			size_t i = 0;
			for (; i < data.size(); ++i)
			{
				const_cast<Eigen::MatrixBase<DI>&>(indices).coeffRef(i) = data[i].index;
				const_cast<Eigen::MatrixBase<DV>&>(values).coeffRef(i) = data[i].value;
			}
			for (; i < nbNeighbours; ++i)
			{
				const_cast<Eigen::MatrixBase<DI>&>(indices).coeffRef(i) = 0;
				const_cast<Eigen::MatrixBase<DV>&>(values).coeffRef(i) = std::numeric_limits<VT>::infinity();
			}
		}
		
#if 0
		//! get the data-point indices from the heap
		/** \return the indices */
		inline IndexVector getIndexes() const
		{
			IndexVector indexes(data.capacity());
			size_t i = 0;
			for (; i < data.size(); ++i)
				indexes.coeffRef(i) = data[i].index;
			for (; i < data.capacity(); ++i)
				indexes.coeffRef(i) = 0;
			return indexes;
		}
#endif
	};
	
#if 0
	//! brute-force implementation of heap
	/** It uses a vector and linear search, which provides replacement in O(n),
	 * 	but with a very low constant overhead. */
	template<typename IT, typename VT>
	struct IndexHeapBruteForceVector
	{
		//! type of an index
		typedef IT Index;
		//! type of a value
		typedef VT Value;
		
		//! an entry of the heap vector
		struct Entry
		{
			IT index; //!< index of point
			VT value;  //!< distance for this point
			
			//! create a new entry
			Entry(const IT index, const VT value): index(index), value(value) {} 
			//! return true if e0 is smaller than e1, false otherwise
			friend bool operator<(const Entry& e0, const Entry& e1) { return e0.value < e1.value; }
		};
		//! vector of entry, type for the storage of the tree
		typedef std::vector<Entry> Entries;
		//! vector of indices
		typedef typename Eigen::Matrix<Index, Eigen::Dynamic, 1> IndexVector;
		
		//! storage for the tree
		Entries data;
		//! reference to the largest value in the tree, to optimise access speed
		const VT& headValueRef;
		//! pre-competed size minus one, to optimise access speed
		const size_t sizeMinusOne;
		
		//! Constructor
		/*! \param size number of elements in the heap */
		IndexHeapBruteForceVector(const size_t size):
			data(size, Entry(0, std::numeric_limits<VT>::infinity())),
			headValueRef((data.end() - 1)->value),
			sizeMinusOne(data.size() - 1)
		{
		}
		
		//! reset to the empty heap
		inline void reset()
		{
			for (typename Entries::iterator it(data.begin()); it != data.end(); ++it)
			{
				it->value = std::numeric_limits<VT>::infinity();
				it->index = 0;
			}
		}
		
		//! get the largest value of the heap
		/** \return the smallest value in the heap */
		inline const VT& headValue() const { return data[0].value; }
		
		//! replace the largest value of the heap
		/** \param index new point index
		 * 	\param value new distance value */
		inline void replaceHead(const Index index, const Value value)
		{
			register size_t i = 0;
			for (; i < sizeMinusOne; ++i)
			{
				if (data[i + 1].value > value)
					data[i] = data[i + 1];
				else
					break;
			}
			data[i].value = value;
			data[i].index = index;
		}
		
		//! sort the entries, from the smallest to the largest
		inline void sort()
		{
			// no need to sort as data are already sorted
		}
		
		//! get the data-point indices from the heap
		/** \return the indices */
		inline IndexVector getIndexes() const
		{
			IndexVector indexes(data.size());
			for (size_t i = 0; i < data.size(); ++i)
				indexes.coeffRef(i) = data[sizeMinusOne-i].index;
			return indexes;
		}
	};
#endif
	
	//! brute-force implementation of heap
	/** It uses a vector and linear search, which provides replacement in O(n),
	 * 	but with a very low constant overhead. */
	template<typename IT, typename VT>
	struct IndexHeapBruteForceVector
	{
		//! type of an index
		typedef IT Index;
		//! type of a value
		typedef VT Value;
		
		//! an entry of the heap vector
		struct Entry
		{
			IT index; //!< index of point
			VT value;  //!< distance for this point
			
			//! create a new entry
			Entry(const IT index, const VT value): index(index), value(value) {} 
			//! return true if e0 is smaller than e1, false otherwise
			friend bool operator<(const Entry& e0, const Entry& e1) { return e0.value < e1.value; }
		};
		//! vector of entry, type for the storage of the tree
		typedef std::vector<Entry> Entries;
		//! vector of indices
		typedef typename Eigen::Matrix<Index, Eigen::Dynamic, 1> IndexVector;
		//! vector of values
		typedef typename Eigen::Matrix<Value, Eigen::Dynamic, 1> ValueVector;
		
		//! storage for the tree
		Entries data;
		//! reference to the largest value in the tree, to optimise access speed
		const VT& headValueRef;
		//! pre-competed size minus one, to optimise access speed
		const size_t sizeMinusOne;
		
		//! Constructor
		/*! \param size number of elements in the heap */
		IndexHeapBruteForceVector(const size_t size):
		data(size, Entry(0, std::numeric_limits<VT>::infinity())),
		headValueRef((data.end() - 1)->value),
		sizeMinusOne(data.size() - 1)
		{
		}
		
		//! reset to the empty heap
		inline void reset()
		{
			for (typename Entries::iterator it(data.begin()); it != data.end(); ++it)
			{
				it->value = std::numeric_limits<VT>::infinity();
				it->index = 0;
			}
		}
		
		//! get the largest value of the heap
		/** \return the smallest value in the heap */
		inline const VT& headValue() const { return headValueRef; }
		
		//! replace the largest value of the heap
		/** \param index new point index
		 * 	\param value new distance value */
		inline void replaceHead(const Index index, const Value value)
		{
			register size_t i;
			for (i = sizeMinusOne; i > 0; --i)
			{
				if (data[i-1].value > value)
					data[i] = data[i-1];
				else
					break;
			}
			data[i].value = value;
			data[i].index = index;
		}
		
		//! sort the entries, from the smallest to the largest
		inline void sort()
		{
			// no need to sort as data are already sorted
		}
		
		//! get the data from the heap
		/** \param indices index vector
		 * 	\param values value vector */
		template<typename DI, typename DV>
		inline void getData(const Eigen::MatrixBase<DI>& indices, const Eigen::MatrixBase<DV> & values) const
		{
			// note: we must implement this hack because of problem with reference to temporary
			// C++0x will solve this with rvalue
			// see: http://eigen.tuxfamily.org/dox-devel/TopicFunctionTakingEigenTypes.html
			// for more informations
			for (size_t i = 0; i < data.size(); ++i)
			{
				const_cast<Eigen::MatrixBase<DI>&>(indices).coeffRef(i) = data[i].index;
				const_cast<Eigen::MatrixBase<DV>&>(values).coeffRef(i) = data[i].value;
			}
		}
#if 0
		//! get the data-point indices from the heap
		/** \return the indices */
		inline IndexVector getIndexes() const
		{
			IndexVector indexes(data.size());
			for (size_t i = 0; i < data.size(); ++i)
				indexes.coeffRef(i) = data[i].index;
			return indexes;
		}
#endif
	};
}

#endif // __INDEX_HEAP_H
