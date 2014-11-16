// Copyright (c) 2014 Romuald PERROT, Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_IMAGE_DESCRIPTION_MLDB_DESCRIPTOR_H
#define OPENMVG_IMAGE_DESCRIPTION_MLDB_DESCRIPTOR_H

#include "openMVG/features/feature.hpp"
#include "openMVG/numeric/math_trait.hpp"

namespace openMVG
{
  /**
    ** @brief Compute final keypoint (ie interest point + description) for a given interest point
    ** @param Li Input Octave slice
    ** @param Lx Input X-derivative
    ** @param Ly Input Y-derivative
    ** @param id_octave Id of current octave
    ** @param ipt Input interest point
    ** @param desc output descriptor (binary descriptor)
    **/
  template< typename Real>
  void ComputeMLDBDescriptor( const Image<Real> & Li , const Image<Real> &Lx , const Image<Real> &Ly ,
                              const int id_octave , const SIOPointFeature & ipt , Descriptor<bool, 486> & desc )
  {
    // TODO
    std::cout << "TODO" << std::endl;
  }
}

#endif
