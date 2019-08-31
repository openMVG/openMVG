// CUDA 加速匹配
#ifndef OPENMVG_MATCHING_CUDA_MATCHER_REGIONS_HPP
#define OPENMVG_MATCHING_CUDA_MATCHER_REGIONS_HPP

#include <memory>

#include "openMVG/matching_image_collection/Matcher.hpp"

namespace openMVG { namespace matching { class PairWiseMatchesContainer; }}
namespace openMVG { namespace sfm { struct Regions_Provider; }}

namespace openMVG {
    namespace matching_image_collection {


        class Cuda_Matcher_Regions : public Matcher {
        public:
            explicit Cuda_Matcher_Regions
                    (
                            float dist_ratio
                    );

            /// Find corresponding points between some pair of view Ids
            void Match
                    (const std::shared_ptr<sfm::Regions_Provider> &regions_provider,
                     const Pair_Set &pairs,
                     matching::PairWiseMatchesContainer &map_PutativesMatches, // the pairwise photometric corresponding points
                     C_Progress *progress = nullptr
                    ) const override;

        private:
            // Distance ratio used to discard spurious correspondence
            float f_dist_ratio_;
        };

    } // namespace matching_image_collection
} // namespace openMVG

#endif // OPENMVG_MATCHING_CUDA_MATCHER_REGIONS_HPP
