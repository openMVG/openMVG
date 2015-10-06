#include "./CasHash.hpp"
#include "openMVG/matching/metric.hpp"
#include "openMVG/stl/indexed_sort.hpp"
#include "nonFree/sift/SIFT_describer.hpp"

// Adapted from author source code:
//
// This source package is aimed to provide a fast and easy-to-use software for
// SIFT keypoints matching in the 3D reconstruction process. For the detailed
// algorithm description, please refer to the CVPR'14 paper.
//
// Fast and Accurate Image Matching with Cascade Hashing for 3D Reconstruction
// Jian Cheng, Cong Leng, Jiaxiang Wu, Hainan Cui, Hanqing Lu. [CVPR 2014]
//
// You are free to use, modify, or redistribute this code in any way you want for
// non-commercial purposes. If you do so, I would appreciate it if you cite the
// CVPR'14 paper.

// Adapted to openMVG in 2014 by cDc and Pierre MOULON under the MPL 2.0 licence.
// Be aware that for commercial purpose you have to ask the CVPR14 authors papers.
// Please cite CVPR'14 paper and openMVG library if you use this implementation.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// D E F I N E S ///////////////////////////////////////////////////

namespace openMVG {
namespace nonFree {
namespace CASHASH {

// S T R U C T S ///////////////////////////////////////////////////

// generate bucket list based on each SIFT point's bucket index
void BuildBuckets(ImageFeatures& imageData)
{
  int cntEleInBucket[kCntBucketPerGroup]; // accumulator; the number of SIFT points in each bucket

  for (int groupIndex = 0; groupIndex < kCntBucketGroup; groupIndex++)
  {
    // initialize <cntEleInBucket>
    memset(cntEleInBucket, 0, sizeof(int)*kCntBucketPerGroup);
    // count the number of SIFT points falling into each bucket
    for (int dataIndex = 0; dataIndex < imageData.cntPoint; dataIndex++)
      cntEleInBucket[imageData.bucketIDList[groupIndex][dataIndex]]++;
    // allocate space for <imageData.bucketList>
    for (int bucketID = 0; bucketID < kCntBucketPerGroup; bucketID++)
    {
      imageData.cntEleInBucket[groupIndex][bucketID] = cntEleInBucket[bucketID];
      imageData.bucketList[groupIndex][bucketID] = new int[cntEleInBucket[bucketID]];
      cntEleInBucket[bucketID] = 0;
    }
    // assign the index of each SIFT point to <imageData.bucketList>
    for (int dataIndex = 0; dataIndex < imageData.cntPoint; dataIndex++)
    {
      const int bucketID = imageData.bucketIDList[groupIndex][dataIndex];
      imageData.bucketList[groupIndex][bucketID][cntEleInBucket[bucketID]++] = dataIndex;
    }
  }
}
/*----------------------------------------------------------------*/


size_t ImportFeatures(
  const std::map<IndexT, std::unique_ptr<features::Regions> > & regions_perImage,
  std::vector<ImageFeatures>& imageDataList)
{
  // accumulator; the number of SIFT feature vectors in <imageDataList>
  size_t cntSiftVec = 0;

  // allocate space for sum and average vectors of SIFT feature
  size_t siftVecSum[kDimSiftData];
  std::vector<int16_t> siftVecAve(kDimSiftData);

  // initialize the sum vector
  memset(siftVecSum, 0, sizeof(size_t)*kDimSiftData);

  // calculate the sum vector by adding up all feature vectors
  for (std::map<IndexT, std::unique_ptr<features::Regions> >::const_iterator iterIma = regions_perImage.begin();
       iterIma != regions_perImage.end(); ++iterIma)
  {
    // Assert we have SIFT regions
    if( dynamic_cast<features::SIFT_Regions*>(regions_perImage.at(0).get()) == NULL )
      return false;

    const DescsT & descriptors = ((features::SIFT_Regions*)regions_perImage.at(0).get())->Descriptors();
    cntSiftVec += descriptors.size();
    for (int dataIndex = 0; dataIndex < descriptors.size(); ++dataIndex)
    {
      const uint8_t* siftVec = descriptors[dataIndex].getData();
      for (int dimIndex = 0; dimIndex < kDimSiftData; ++dimIndex)
      {
        siftVecSum[dimIndex] += *(siftVec++);
      }
    }
  }

  // calculate the average vector
  for (int dimIndex = 0; dimIndex < kDimSiftData; dimIndex++)
    siftVecAve[dimIndex] = (int16_t)(siftVecSum[dimIndex] / cntSiftVec);

  // init feature descriptors for each image
  imageDataList.resize(regions_perImage.size());
  const HashConvertor stHashConvertor;
#ifdef OPENMVG_USE_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (int i = 0; i < regions_perImage.size(); ++i)
  {
    std::map<IndexT, std::unique_ptr<features::Regions> >::const_iterator iterDesc = regions_perImage.begin();
    std::advance(iterDesc, i);

    const DescsT & descriptors = ((features::SIFT_Regions*)iterDesc->second.get())->Descriptors();

    ImageFeatures& imageDataDst = imageDataList[i];
    // substract the average vector from each feature vector
    imageDataDst.cntPoint = descriptors.size();

    // convert SIFT feature to Hash code
    stHashConvertor.SiftDataToHashData(imageDataDst, descriptors, siftVecAve);

    // construct multiple groups of buckets for each image
    BuildBuckets(imageDataDst);
  }

  return cntSiftVec;
}
/*----------------------------------------------------------------*/

// C L A S S  //////////////////////////////////////////////////////

HashConvertor::HashConvertor()
{
  // generate the projection matrix for the primary hashing function (for Hamming distance ranking)
  for (int index_1 = 0; index_1 < kDimHashData; ++index_1)
  {
    for (int index_2 = 0; index_2 < kDimSiftData; ++index_2)
    {
      projMatPriTr[index_1][index_2] = static_cast<int>(GetNormRand() * 1000);
    }
  }

#ifdef Bucket_SecHash
  // generate the projection matrix for the secondary hashing function (for bucket construction)
  for (int groupIndex = 0; groupIndex < kCntBucketGroup; ++groupIndex)
  {
    for (int index_1 = 0; index_1 < kCntBucketBit; ++index_1)
    {
      for (int index_2 = 0; index_2 < kDimSiftData; ++index_2)
      {
        projMatSecTr[groupIndex][index_1][index_2] = static_cast<int>(GetNormRand() * 1000);
      }
    }
  }
#endif // Bucket_SecHash

#ifdef Bucket_PriHashSel
  // generate the selected bit list (for bucket construction)
  // bits are selected from the Hash code generated by the primary hashing function
  for (int groupIndex = 0; groupIndex < kCntBucketGroup; groupIndex++)
  {
    uint8_t bitUsedList[kDimHashData];
    memset(bitUsedList, 0, kDimHashData); // initialize the bit usage flag array

    for (int bitIndex = 0; bitIndex < kCntBucketBit; bitIndex++)
    {
      int bitSel = -1;
      do
      {
        bitSel = rand() % kDimHashData;
      }   while (bitUsedList[bitSel] == 1); // ensure the selected bit has not been used for this bucket

      bitUsedList[bitSel] = 1;
      bucketBitList[groupIndex][bitIndex] = bitSel;
    }
  }
#endif // Bucket_PriHashSel
}

void HashConvertor::SiftDataToHashData(
  ImageFeatures& imageData,
  const DescsT & descs,
  std::vector<int16_t> & averageDescriptor) const
{
  // allocate space for <compHashDataPtrList> and <bucketIDList>
  std::vector<HashData> hashDataList(imageData.cntPoint); // Hash code for each SIFT point
  imageData.compHashDataPtrList = new CompHashDataPtr[imageData.cntPoint];
  for (int groupIndex = 0; groupIndex < kCntBucketGroup; groupIndex++)
    imageData.bucketIDList[groupIndex] = new uint16_t[imageData.cntPoint];

  std::vector<int16_t> tempDesc(kDimSiftData);
  for (int dataIndex = 0; dataIndex < imageData.cntPoint; ++dataIndex)
  {
    // allocate space for each SIFT point
    imageData.compHashDataPtrList[dataIndex] = new uint64_t[kDimCompHashData];

    // obtain pointers for SIFT feature vector, Hash code and CompHash code

    // -- compute on the fly the mean centered SIFT feature vector
    const DescriptorT & desc = descs[dataIndex];
    for (int dimIndex = 0; dimIndex < kDimSiftData; ++dimIndex) {
      tempDesc[dimIndex] = desc[dimIndex] - averageDescriptor[dimIndex];
    }

    const SiftData* siftDataPtr = &tempDesc[0];
    HashData& hashData = hashDataList[dataIndex];
    CompHashDataPtr compHashDataPtr = imageData.compHashDataPtrList[dataIndex];

    // calculate the Hash code: H = PX
    // H: Hash code vector, Nx1 vector, N = kDimHashData
    // P: projection matrix, NxM matrix, M = kDimSiftData
    // X: SIFT feature vector, Mx1 vector
    for (int dimHashIndex = 0; dimHashIndex < kDimHashData; dimHashIndex++)
    {
      int sum = 0;
      const int* projVec = projMatPriTr[dimHashIndex];
      for (int dimSiftIndex = 0; dimSiftIndex < kDimSiftData; dimSiftIndex++)
      {
        if (siftDataPtr[dimSiftIndex] != 0)
        {
          sum += (int)siftDataPtr[dimSiftIndex] * projVec[dimSiftIndex]; // matrix multiplication
        }
      }
      hashData.hash[dimHashIndex] = (sum > 0) ? 1 : 0; // use 0 as the threshold
    }

    // calculate the CompHash code
    // compress <kBitInCompHash> Hash code bits within a single <uint64_t> variable
    for (int dimCompHashIndex = 0; dimCompHashIndex < kDimCompHashData; dimCompHashIndex++)
    {
      uint64_t compHashBitVal = 0;
      const int dimHashIndexLBound = dimCompHashIndex * kBitInCompHash;
      const int dimHashIndexUBound = (dimCompHashIndex + 1) * kBitInCompHash;
      for (int dimHashIndex = dimHashIndexLBound; dimHashIndex < dimHashIndexUBound; ++dimHashIndex)
        compHashBitVal = (compHashBitVal << 1) + hashData.hash[dimHashIndex]; // set the corresponding bit to 1/0
      compHashDataPtr[dimCompHashIndex] = compHashBitVal;
    }

    // determine the bucket index for each bucket group
    for (int groupIndex = 0; groupIndex < kCntBucketGroup; ++groupIndex)
    {
      uint16_t bucketID = 0;
      for (int bitIndex = 0; bitIndex < kCntBucketBit; ++bitIndex)
      {
#ifdef Bucket_SecHash
        int sum = 0;
        const int* projVec = projMatSecTr[groupIndex][bitIndex];
        for (int dimSiftIndex = 0; dimSiftIndex < kDimSiftData; dimSiftIndex++)
        {
          if (siftDataPtr[dimSiftIndex] != 0)
          {
            sum += (int)siftDataPtr[dimSiftIndex] * projVec[dimSiftIndex]; // matrix multiplication
          }
        }
        bucketID = (bucketID << 1) + (sum > 0 ? 1 : 0);
#endif // Bucket_SecHash
#ifdef Bucket_PriHashSel
        bucketID = (bucketID << 1) + hashData.hash[bucketBitList[groupIndex][bitIndex]];
#endif // Bucket_PriHashSel
      }
      imageData.bucketIDList[groupIndex][dataIndex] = bucketID;
    }
  }
}

double HashConvertor::GetNormRand()
{
  // based on Box-Muller transform; for more details, please refer to the following WIKIPEDIA website:
  // http://en.wikipedia.org/wiki/Box_Muller_transform
  const double u1 = (rand() % 1000 + 1) / 1000.0;
  const double u2 = (rand() % 1000 + 1) / 1000.0;
  return sqrt(-2.0 * log(u1)) * cos(2.0 * acos(-1.0) * u2);
}
/*----------------------------------------------------------------*/



// C L A S S  //////////////////////////////////////////////////////

template<>
void CasHashMatcher::MatchSpFast<DescsT>(
  std::vector<matching::IndMatch>& matchList,
  const ImageFeatures& imageData_1,
  const DescsT & descriptors1,
  const ImageFeatures& imageData_2,
  const DescsT & descriptors2,
  float ratiomax) const
{
  typedef openMVG::matching::L2_Vectorized<DescriptorT::bin_type> MetricL2T;
  MetricL2T MetricL2;

  std::vector<int> candidateIndexList;  // indexes of candidate SIFT points
  std::set<int> dataIndexUsedList;      // usage flag list of SIFT points; to indicate whether this point has been added to <candidateIndexListTop>
  std::vector<unsigned char> distList;  // Hamming distance of candidate SIFT points to the query point
  std::vector<int> linkList[kDimHashData+1]; // re-assign candidate SIFT points according to their Hamming distance (the number of candidates with Hamming distance of [0, 1, ..., 128])

  // reserve some memory for the temporary arrays
  const size_t maxCount = std::max(descriptors1.size(), descriptors2.size());
  candidateIndexList.reserve(maxCount);
  distList.reserve(maxCount);

  int candidateIndexListTop[kCntCandidateTopMax]; // indexes of final candidate SIFT points
  MetricL2T::ResultType candidateDistListTop[kCntCandidateTopMax]; // Euclidean distance of final candidate SIFT points

  const float ratiomaxSq(Square(ratiomax));

  // initialize <matchList>
  matchList.clear();

  // try to find the corresponding SIFT point in <imageData_2> for each point in <imageData_1>
  for (int dataIndex_1 = 0; dataIndex_1 < imageData_1.cntPoint; ++dataIndex_1)
  {
    candidateIndexList.clear();
    dataIndexUsedList.clear();

    // fetch candidate SIFT points from the buckets in each group
    // only the bucket with the same <bucketID> with the current query point will be fetched

    for (int groupIndex = 0; groupIndex < kCntBucketGroup; ++groupIndex)
    {
      // obtain the <bucketID> for the query SIFT point
      const int bucketID = imageData_1.bucketIDList[groupIndex][dataIndex_1];
      // obtain the pointer to the corresponding bucket
      const int* bucketPtr = imageData_2.bucketList[groupIndex][bucketID];

      for (int eleIndex = 0; eleIndex < imageData_2.cntEleInBucket[groupIndex][bucketID]; ++eleIndex)
      {
        candidateIndexList.push_back(bucketPtr[eleIndex]); // fetch candidates from the bucket
      }
    }
    const size_t cntCandidate = candidateIndexList.size();
    distList.resize(cntCandidate);
    // calculate the Hamming distance of all candidates based on the CompHash code
    const CompHashDataPtr ptr_1 = imageData_1.compHashDataPtrList[dataIndex_1];
    for (int candidateIndex = 0; candidateIndex < cntCandidate; ++candidateIndex)
    {
      const CompHashDataPtr ptr_2 = imageData_2.compHashDataPtrList[candidateIndexList[candidateIndex]];
      distList[candidateIndex] =
        openMVG::matching::Hamming<uint64_t>::popcnt64(ptr_1[0] ^ ptr_2[0]) +
        openMVG::matching::Hamming<uint64_t>::popcnt64(ptr_1[1] ^ ptr_2[1]);
    }

    // re-assign candidates to a linked list based on their Hamming distance
    for (int i = 0; i <= kDimHashData; ++i)
      linkList[i].clear();
    for (int candidateIndex = 0; candidateIndex < cntCandidate; ++candidateIndex)
    {
      const uint8_t hashDist = distList[candidateIndex];
      linkList[hashDist].push_back(candidateIndexList[candidateIndex]);
    }

    // add top-ranked candidates in Hamming distance to <candidateIndexListTop>
    int cntCandidateFound = 0;
    for (int hashDist = 0;
         hashDist <= kDimHashData &&
         cntCandidateFound < kCntCandidateTopMin; // enough candidates have been selected: break
         ++hashDist)
    {
      for (int linkListIndex = linkList[hashDist].size()-1;
           linkListIndex >= 0 &&
           cntCandidateFound < kCntCandidateTopMax; //enough candidates have been selected: break
           --linkListIndex)
      {
        const int dataIndex_2 = linkList[hashDist][linkListIndex];
        if (dataIndexUsedList.find(dataIndex_2) == dataIndexUsedList.end()) // avoid selecting same candidate multiple times
        {
          dataIndexUsedList.insert(dataIndex_2);
          candidateIndexListTop[cntCandidateFound++] = dataIndex_2; // add candidate to <candidateIndexListTop>
        }
      }
    }

    // calculate Euclidean distance for candidates in <candidateIndexListTop>
    for (int candidateIndex = 0; candidateIndex < cntCandidateFound; ++candidateIndex)
    {
      const int dataIndex_2 = candidateIndexListTop[candidateIndex];

      // fetch the pointers to the two chosen SIFT feature vectors
      const DescsT::value_type::bin_type * ptr_1 = descriptors1[dataIndex_1].getData();
      const DescsT::value_type::bin_type * ptr_2 = descriptors2[dataIndex_2].getData();
      candidateDistListTop[candidateIndex] = MetricL2(ptr_1, ptr_2, DescsT::value_type::static_size);
    }

    // find the top-2 candidates with minimal Euclidean distance
    MetricL2T::ResultType minVal_1(0), minVal_2(0);
    int minValInd_1 = -1, minValInd_2 = -1;

    for (int candidateIndex = 0; candidateIndex < cntCandidateFound; ++candidateIndex)
    {
      if (minValInd_2 == -1 || minVal_2 > candidateDistListTop[candidateIndex])
      {
        minVal_2 = candidateDistListTop[candidateIndex];
        minValInd_2 = candidateIndexListTop[candidateIndex];
      }
      if (minValInd_1 == -1 || minVal_1 > minVal_2)
      {
        std::swap(minVal_1, minVal_2);
        std::swap(minValInd_1, minValInd_2);
      }
    }

    // apply the threshold for matching rejection
    if (cntCandidateFound >=2 && minVal_1 < ratiomaxSq * minVal_2)
      matchList.push_back(openMVG::matching::IndMatch(dataIndex_1, minValInd_1));
  }
  if (matchList.size() < kMinMatchListLen)
    matchList.clear();
}
/*----------------------------------------------------------------*/

} // namespace CASHASH
} // namespace nonFree
} // namespace openMVG
