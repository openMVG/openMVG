//
// Created by saurus on 18-7-26.
//

#ifndef OPENMVG_SFM_DATA_PBA_H
#define OPENMVG_SFM_DATA_PBA_H

#include <unordered_map>
#include <third_party/pba/src/pba/pba.h>
#include <openMVG/cameras/Camera_Common.hpp>
#include "sfm_data.hpp"

namespace openMVG {
  namespace sfm {
    class sfm_data_PBA {
    public:
      sfm_data_PBA(
        ParallelBA::DeviceT device
      ) {
      }

      bool DataToPBA(SfM_Data &sfm_data,
                     cameras::Intrinsic_Parameter_Type intrinsic_refinement_options);
      bool Adjust(SfM_Data &sfm_data);

    private:
      // pba input & output
      vector<CameraT> camera_data;
      vector<Point3D> point_data;

      // pba input
      vector<Point2D> measurements;
      vector<int> camidx, ptidx;
      vector<int> focalmask;

      ParallelBA pba;
      std::unordered_map<unsigned long, unsigned long> view_id2camera_id;
    };
  }
}


#endif //OPENMVG_SFM_DATA_PBA_H
