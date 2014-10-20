#ifndef _CASHASH_H_
#define _CASHASH_H_

// Adapted from author source code:
//
// This source package is aimed to provide a fast and easy-to-use software for
// SIFT keypoints matching in the 3D reconstruction process. For the detailed
// algorithm description, please refer to our CVPR'14 paper.
//
// Fast and Accurate Image Matching with Cascade Hashing for 3D Reconstruction
// Jian Cheng, Cong Leng, Jiaxiang Wu, Hainan Cui, Hanqing Lu. [CVPR 2014]
//
// You are free to use, modify, or redistribute this code in any way you want for
// non-commercial purposes. If you do so, I would appreciate it if you cite our
// CVPR'14 paper.

// Adapted to openMVG in 2014 by cDc and Pierre MOULON under the MPL 2.0 licence.
// Be aware that for commercial purpose you have to ask the CVPR14 authors papers.
// Please cite CVPR'14 paper and openMVG library if you use this implementation.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// I N C L U D E S /////////////////////////////////////////////////

#include "openMVG/features/features.hpp"
#include "openMVG/matching/indMatch.hpp"
#include <stdint.h>

// D E F I N E S ///////////////////////////////////////////////////

#define Bucket_SecHash // use secondary hashing function to construct buckets
//#define Bucket_PriHashSel // use selected bits in primary hashing function to construct buckets

// P R O T O T Y P E S /////////////////////////////////////////////

namespace openMVG {
namespace nonFree {
namespace CASHASH {

const int kDimSiftData = 128; // the number of dimensions of SIFT feature
const int kDimHashData = 128; // the number of dimensions of Hash code
const int kBitInCompHash = 64; // the number of Hash code bits to be compressed; in this case, use a <uint64_t> variable to represent 64 bits
const int kDimCompHashData = kDimHashData / kBitInCompHash; // the number of dimensions of CompHash code
const int kMinMatchListLen = 16; // the minimal list length for outputing SIFT matching result between two images

const int kCntBucketBit = 8; //10; // the number of bucket bits
const int kCntBucketGroup = 6; // the number of bucket groups
const int kCntBucketPerGroup = 1 << kCntBucketBit; // the number of buckets in each group

const int kCntCandidateTopMin = 6; // the minimal number of top-ranked candidates
const int kCntCandidateTopMax = 10; // the maximal number of top-ranked candidates

typedef int16_t SiftData; // SIFT feature components are represented with <int16_t> type
struct HashData {
  uint8_t hash[kDimHashData];
}; // Hash code is represented with <uint8_t> type; only the lowest bit is used
typedef uint64_t* CompHashDataPtr; // CompHash code is represented with <uint64_t> type
typedef int* BucketElePtr; // index list of points in a specific bucket

// all information needed for an image to perform CasHash-Matching
struct ImageFeatures {
  int cntPoint; // the number of SIFT points
  CompHashDataPtr* compHashDataPtrList; // CompHash code for each SIFT point
  uint16_t* bucketIDList[kCntBucketGroup]; // bucket entries for each SIFT point
  int cntEleInBucket[kCntBucketGroup][kCntBucketPerGroup]; // the number of SIFT points in each bucket
  BucketElePtr bucketList[kCntBucketGroup][kCntBucketPerGroup]; // SIFT point index list for all buckets

  ImageFeatures() {
    memset(this, 0, sizeof(ImageFeatures));
  }
  ~ImageFeatures() {
    for (int groupIndex = 0; groupIndex < kCntBucketGroup; ++groupIndex) {
      delete[] bucketIDList[groupIndex];
      for (int bucketID = 0; bucketID < kCntBucketPerGroup; ++bucketID)
        delete[] bucketList[groupIndex][bucketID];
    }
    for (int dataIndex = 0; dataIndex < cntPoint; ++dataIndex)
      delete[] compHashDataPtrList[dataIndex];
    delete[] compHashDataPtrList;
  }
};
/*----------------------------------------------------------------*/

// Define SIFT descriptors
typedef Descriptor<unsigned char, 128> DescriptorT;
typedef std::vector<DescriptorT > DescsT;
typedef std::map<size_t, DescsT > map_DescT;

// fetch and adjust each SIFT feature vector in <imageDataList> to zero-mean
// returns the total number of features
size_t ImportFeatures(const map_DescT& map_Desc, std::vector<ImageFeatures>& imageDataList);
/*----------------------------------------------------------------*/


// this class is used to convert SIFT feature vectors to Hash codes and CompHash codes
class HashConvertor
{
public:
  // constructor function; to initialize private data member variables
  HashConvertor();
  // convert SIFT feature vectors to Hash codes and CompHash codes
  // also, the bucket index for each SIFT point will be determined (different bucket groups correspond to different indexes)
  void SiftDataToHashData(
    ImageFeatures& imageData,
    const DescsT & descs,
    std::vector<int16_t> & averageDescriptor) const;

private:
  int projMatPriTr[kDimHashData][kDimSiftData]; // projection matrix of the primary hashing function
#ifdef Bucket_SecHash
  int projMatSecTr[kCntBucketGroup][kCntBucketBit][kDimSiftData]; // projection matrix of the secondary hashing function
#endif // Bucket_SecHash
#ifdef Bucket_PriHashSel
  int bucketBitList[kCntBucketGroup][kCntBucketBit]; // selected bits in the result of primary hashing function for bucket construction
#endif // Bucket_PriHashSel

private:
  // generate random number which follows normal distribution, with <mean = 0> and <variance = 1>
  static double GetNormRand();
};
/*----------------------------------------------------------------*/


// this class is used to perform SIFT points matching using CasHash Matching Technique
class CasHashMatcher
{
public:
  // return a list of matched SIFT points between two input images (ratiomax - maximum distance ratio)
  template<typename Desc>
  void MatchSpFast(
    std::vector<matching::IndMatch>& matchList,
    const ImageFeatures& imageData_1,
    const Desc & descriptors1,
    const ImageFeatures& imageData_2,
    const Desc & descriptors2,
    float ratiomax = 0.8f) const
  {
    std::cout << "Not implemented for this Descriptor type" << std::endl;
  }
};

template<>
void CasHashMatcher::MatchSpFast<DescsT>(
  std::vector<matching::IndMatch>& matchList,
  const ImageFeatures& imageData_1,
  const DescsT & descriptors1,
  const ImageFeatures& imageData_2,
  const DescsT & descriptors2,
  float ratiomax) const;
/*----------------------------------------------------------------*/

} // namespace CASHASH
} // namespace nonFree
} // namespace openMVG

#endif //_CASHASH_H_
