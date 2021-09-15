// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_PAIRWISE_ADJACENCY_DISPLAY_HPP
#define OPENMVG_PAIRWISE_ADJACENCY_DISPLAY_HPP

#include <string>

#include "openMVG/matching/indMatch.hpp"
#include "openMVG/graphics/color_gradient.hpp"
#include "third_party/vectorGraphics/svgDrawer.hpp"

namespace openMVG  {
namespace matching {

/// Display pair wises matches as an Adjacency matrix in svg format
void PairWiseMatchingToAdjacencyMatrixSVG
(
  const size_t NbImages,
  const matching::PairWiseMatches & map_Matches,
  const std::string & sOutName
)
{
  if ( !map_Matches.empty())
  {
    // Set the coloring gradient interface
    graphics::Color_Gradient heatMapGradient(graphics::Color_Gradient::k2BlueRedHeatMap());
    float max_match_count = 0;
    for (const auto & match_it : map_Matches)
    {
      max_match_count = std::max(max_match_count, static_cast<float>(match_it.second.size()));
    }

    const float scaleFactor = 5.0f;
    svg::svgDrawer svgStream((NbImages+3)*5, (NbImages+3)*5);
    // Go along all possible pair
    for (size_t I = 0; I < NbImages; ++I) {
      for (size_t J = 0; J < NbImages; ++J) {
        // If the pair have matches display a blue boxes at I,J position.
        auto iterSearch = map_Matches.find({I,J});
        if (iterSearch != map_Matches.end() && !iterSearch->second.empty())
        {
          // Display as a tooltip: "(IndexI, IndexJ NbMatches)"
          std::ostringstream os_tooltip;
          os_tooltip << "(" << J << "," << I << " " << iterSearch->second.size() <<")";

          float r,g,b;
          heatMapGradient.getColor(iterSearch->second.size() / max_match_count, r, g, b);
          std::ostringstream os_color;
          os_color << "rgb(" << int(r * 255) << "," << int(g  * 255) << "," << int(b * 255) << ")";

          svgStream.drawSquare(J*scaleFactor, I*scaleFactor, scaleFactor/2.0f,
            svg::svgStyle().fill(os_color.str()).noStroke().tooltip(os_tooltip.str()));
        }
      }
    }
    // Display axes with 0 -> NbImages annotation : _|
    std::ostringstream osNbImages;
    osNbImages << NbImages;
    svgStream.drawText((NbImages+1)*scaleFactor, scaleFactor, scaleFactor, "0", "black");
    svgStream.drawText((NbImages+1)*scaleFactor,
      (NbImages)*scaleFactor - scaleFactor, scaleFactor, osNbImages.str(), "black");
    svgStream.drawLine((NbImages+1)*scaleFactor, 2*scaleFactor,
      (NbImages+1)*scaleFactor, (NbImages)*scaleFactor - 2*scaleFactor,
      svg::svgStyle().stroke("black", 1.0));

    svgStream.drawText(scaleFactor, (NbImages+1)*scaleFactor, scaleFactor, "0", "black");
    svgStream.drawText((NbImages)*scaleFactor - scaleFactor,
      (NbImages+1)*scaleFactor, scaleFactor, osNbImages.str(), "black");
    svgStream.drawLine(2*scaleFactor, (NbImages+1)*scaleFactor,
      (NbImages)*scaleFactor - 2*scaleFactor, (NbImages+1)*scaleFactor,
      svg::svgStyle().stroke("black", 1.0));

    std::ofstream svgFileStream( sOutName.c_str());
    svgFileStream << svgStream.closeSvgFile().str();
  }
}

} // namespace matching
} // namespace openMVG

#endif // OPENMVG_PAIRWISE_ADJACENCY_DISPLAY_HPP
