// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/matching_image_collection/Cuda_Matcher_Regions.hpp"

#include "openMVG/matching/cascade_hasher.hpp"
#include "openMVG/features/feature.hpp"
#include "openMVG/features/image_describer.hpp"
#include "openMVG/features/regions_factory.hpp"
#include "openMVG/matching/matching_filters.hpp"
#include "openMVG/matching/indMatchDecoratorXY.hpp"
#include "openMVG/sfm/pipelines/sfm_regions_provider.hpp"
#include "openMVG/types.hpp"
#include <Eigen/Core>

#include "third_party/cudasift/cudaSift.h"

#include "third_party/progress/progress.hpp"


namespace openMVG {
    namespace matching_image_collection {

        using namespace openMVG::matching;
        using namespace openMVG::features;

        void FillSiftData(SiftData &data, const std::shared_ptr<features::Regions> regions) {
            const std::vector<features::PointFeature> pointFeatures = regions->GetRegionsPositions();
            /*const ScalarT *tab =
                    reinterpret_cast<const ScalarT *>(regions->DescriptorRawData());*/
            const size_t dimension = regions->DescriptorLength();
            int pt_count = pointFeatures.size();

            SIFT_Regions *regionsCasted = dynamic_cast<SIFT_Regions *>(regions.get());

            size_t region_count = regionsCasted->RegionCount();

            InitSiftData(data, region_count, true, true);
            data.numPts = region_count;
#ifdef MANAGEDMEM
            SiftPoint *h_data = data.m_data;
#else
            SiftPoint *h_data = data.h_data;
#endif

#ifdef OPENMVG_USE_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
            for (int i = 0; i < region_count; i++) {
                const SIOPointFeature &feature = regionsCasted->Features()[i];
                Descriptor<unsigned char, 128> descriptor = regionsCasted->Descriptors()[i];

                h_data[i].xpos = feature.x();
                h_data[i].ypos = feature.y();
                h_data[i].scale = feature.scale();
                h_data[i].orientation = feature.orientation();

                for (unsigned int ctr = 0; ctr < 128; ++ctr) {
                    h_data[i].data[ctr] = static_cast<float>(descriptor[ctr]) / 512.0f;
                }
            }
            CopySiftDataHost2Dev(data);
        }


        void UpdateMatch(SiftData &siftData1, SiftData &siftData2, float *homography)
        {
#ifdef MANAGEDMEM
            SiftPoint *sift1 = siftData1.m_data;
  SiftPoint *sift2 = siftData2.m_data;
#else
            SiftPoint *sift1 = siftData1.h_data;
            SiftPoint *sift2 = siftData2.h_data;
#endif
            int numPts1 = siftData1.numPts;
            int numPts2 = siftData2.numPts;
            int numFound = 0;
#if 0
            homography[0] = homography[4] = -1.0f;
            homography[1] = homography[3] = homography[6] = homography[7] = 0.0f;
            homography[2] = 1279.0f;
            homography[5] = 959.0f;
#endif

#ifdef OPENMVG_USE_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
            for (int i=0;i<numPts1;i++) {
                float *data1 = sift1[i].data;
                float den = homography[6]*sift1[i].xpos + homography[7]*sift1[i].ypos + homography[8];
                float dx = (homography[0]*sift1[i].xpos + homography[1]*sift1[i].ypos + homography[2]) / den - sift1[i].match_xpos;
                float dy = (homography[3]*sift1[i].xpos + homography[4]*sift1[i].ypos + homography[5]) / den - sift1[i].match_ypos;
                float err = dx*dx + dy*dy;
                sift1[i].match_error = err;
            }

            return;

            for (int i=0;i<numPts1;i++) {
                float *data1 = sift1[i].data;
                //std::cout << i << ":" << sift1[i].scale << ":" << (int)sift1[i].orientation << " " << sift1[i].xpos << " " << sift1[i].ypos << std::endl;
                bool found = false;
                for (int j=0;j<numPts2;j++) {
                    if(sift1[i].match == j)
                    {
                        continue;
                    }

                    float *data2 = sift2[j].data;
                    float sum = 0.0f;
                   /* for (int k=0;k<128;k++) {
                        sum += data1[k] * data2[k];
                    }*/
                    float den = homography[6]*sift1[i].xpos + homography[7]*sift1[i].ypos + homography[8];
                    float dx = (homography[0]*sift1[i].xpos + homography[1]*sift1[i].ypos + homography[2]) / den - sift2[j].xpos;
                    float dy = (homography[3]*sift1[i].xpos + homography[4]*sift1[i].ypos + homography[5]) / den - sift2[j].ypos;
                    float err = dx*dx + dy*dy;
                    found = err<100.0f;
                    if(found && err < sift1[i].match_error)
                    {
                        sift1[i].match = j;
                        sift1[i].match_xpos = sift2[j].xpos;
                        sift1[i].match_ypos = sift2[j].ypos;
                        sift1[i].match_error = err;
                    }
                }
            }
        }

        Cuda_Matcher_Regions::Cuda_Matcher_Regions(float distRatio)
                : Matcher(), f_dist_ratio_(distRatio) {
        }

        namespace impl2 {
            template<typename ScalarT>
            void Match(const sfm::Regions_Provider &regions_provider,
                       const Pair_Set &pairs,
                       float fDistRatio,
                       PairWiseMatchesContainer &map_PutativesMatches, // the pairwise photometric corresponding points
                       C_Progress *my_progress_bar) {
                if (!my_progress_bar)
                    my_progress_bar = &C_Progress::dummy();
                my_progress_bar->restart(pairs.size(), "\n- Matching -\n");

                // Collect used view indexes
                std::set<IndexT> used_index;
                // Sort pairs according the first index to minimize later memory swapping
                using Map_vectorT = std::map<IndexT, std::vector<IndexT>>;
                Map_vectorT map_Pairs;
                for (const auto &pair_idx : pairs) {
                    map_Pairs[pair_idx.first].push_back(pair_idx.second);
                    used_index.insert(pair_idx.first);
                    used_index.insert(pair_idx.second);
                }

                SiftData siftData_I, siftData_J;

                // Perform matching between all the pairs
                for (const auto &pairs : map_Pairs) {
                    if (my_progress_bar->hasBeenCanceled())
                        break;
                    const IndexT I = pairs.first;
                    const std::vector<IndexT> &indexToCompare = pairs.second;

                    const std::shared_ptr<features::Regions> regionsI = regions_provider.get(I);
                    if (regionsI->RegionCount() == 0) {
                        (*my_progress_bar) += indexToCompare.size();
                        continue;
                    }

                    const std::vector<features::PointFeature> pointFeaturesI = regionsI->GetRegionsPositions();

                    // 开始构造SiftData

                    FillSiftData(siftData_I, regionsI);
                    //std::cout << "构建SiftData I完成; 特征数量: "<< siftData_I.numPts << " ;" << std::endl;


                    /*#ifdef OPENMVG_USE_OPENMP
                    #pragma omp parallel for schedule(dynamic)
                    #endif*/
                    for (int j = 0; j < (int)indexToCompare.size(); ++j) {
                        if (my_progress_bar->hasBeenCanceled())
                            continue;
                        const size_t J = indexToCompare[j];
                        const std::shared_ptr<features::Regions> regionsJ = regions_provider.get(J);

                        // 开始构造SiftData
                        FillSiftData(siftData_J, regionsJ);
                        //std::cout << "构建SiftData J完成; 特征数量: "<< siftData_I.numPts << " ;" << std::endl;

                        double match_time = MatchSiftData(siftData_I, siftData_J);

                        float homography[9];
                        int numMatches;
                        FindHomography(siftData_I, homography, &numMatches, 10000, 0.00f, 0.90f, 3.5);

                        /* 打印匹配信息
                        std::cout << "Number of original features: " <<  siftData_I.numPts << " " << siftData_J.numPts << std::endl;
                        std::cout << "Number of matching features: "<< numMatches << " "
                        << 100.0f * numMatches/std::min(siftData_I.numPts, siftData_J.numPts) << std::endl;
                         */

                        UpdateMatch(siftData_I,siftData_J,homography);

                        matching::IndMatches vec_putative_matches;

                        float  sum_err = 0.0;
                        int num_mt=0;
/*#ifdef OPENMVG_USE_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif*/
                        for(int m=0; m< siftData_I.numPts; ++m)
                        {
                            if(siftData_I.h_data[m].match_error <100) {
                                ++num_mt;
                                int match = siftData_I.h_data[m].match;
                                float match_err = siftData_I.h_data[m].match_error;
                                sum_err += match_err;

                                vec_putative_matches.emplace_back( m , match);


                                /* 输出匹配点
                                 * std::cout << m << ": (" << siftData_I.h_data[m].xpos
                                          << ", " << siftData_I.h_data[m].ypos
                                          << "), "<< match <<": (" << siftData_I.h_data[m].match_xpos
                                          << ", " << siftData_I.h_data[m].match_ypos
                                          << "), Error: " << match_err
                                          << std::endl;
                                */
                            }
                        }

                        std::cout << "Pair: "<<I <<" - "<<J<< ", Number of matching features: "<< num_mt <<"; ";
                        /* 统计误差 */
                        if(num_mt>0) {
                            float rmse = sqrt(sum_err / num_mt);
                            std::cout << "RMSE: " << rmse;
                        }
                        std::cout << std::endl;

//                        std::cout << "匹配 "<< I <<" - "<< J <<" 耗时: "<< match_time << " ;" << std::endl;
//                        std::cout << "  "<< I <<" 匹配点数: "<< i_match_num <<" 误差: "<< i_match_err << " ;" << std::endl;
//                        std::cout << "  "<< J <<" 匹配点数: "<< j_match_num <<" 耗时: "<< j_match_err << " ;" << std::endl;

                        FreeSiftData(siftData_J);

#ifdef OPENMVG_USE_OPENMP
#pragma omp critical
#endif
                        {
                            if (!vec_putative_matches.empty())
                            {
                                map_PutativesMatches.insert(
                                        {
                                                {I,J},
                                                std::move(vec_putative_matches)
                                        });
                            }
                        }
                        ++(*my_progress_bar);

                    }
                    FreeSiftData(siftData_I);
                }
            }
        } // namespace impl

        void Cuda_Matcher_Regions::Match
                (
                        const std::shared_ptr<sfm::Regions_Provider> &regions_provider,
                        const Pair_Set &pairs,
                        PairWiseMatchesContainer &map_PutativesMatches, // the pairwise photometric corresponding points
                        C_Progress *my_progress_bar
                ) const {
#ifdef OPENMVG_USE_OPENMP
            std::cout << "Using the OPENMP thread interface" << std::endl;
#endif
            if (!regions_provider)
                return;

            if (regions_provider->IsBinary())
                return;

            if (regions_provider->Type_id() == typeid(unsigned char).name()) {
                impl2::Match<unsigned char>(
                        *regions_provider.get(),
                        pairs,
                        f_dist_ratio_,
                        map_PutativesMatches,
                        my_progress_bar);
            } else if (regions_provider->Type_id() == typeid(float).name()) {
                impl2::Match<float>(
                        *regions_provider.get(),
                        pairs,
                        f_dist_ratio_,
                        map_PutativesMatches,
                        my_progress_bar);
            } else {
                std::cerr << "Matcher not implemented for this region type" << std::endl;
            }
        }

    } // namespace openMVG
} // namespace matching_image_collection
