//
// Created by saurus on 18-7-26.
//

#ifndef OPENMVG_SFM_DATA_PBA_H
#define OPENMVG_SFM_DATA_PBA_H

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
      vector<CameraT> camera_data;
      vector<Point2D> measurements;
      vector<Point3D> point_data;
      vector<int> camidx, ptidx;
      ParallelBA pba;
      int *focalmask;
    };
  }
}


#endif //OPENMVG_SFM_DATA_PBA_H
