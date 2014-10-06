
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_PAIRWISE_ADJACENCY_DISPLAY_H
#define OPENMVG_PAIRWISE_ADJACENCY_DISPLAY_H

#include <algorithm>

#include "third_party/vectorGraphics/svgDrawer.hpp"
using namespace svg;

#include "openMVG/matching/indMatch.hpp"
using namespace openMVG::matching;

namespace openMVG  {

/// Display pair wises matches as an Adjacency matrix in svg format
void PairWiseMatchingToAdjacencyMatrixSVG(const size_t NbImages,
  const matching::PairWiseMatches & map_Matches,
  const std::string & sOutName)
{
  if ( !map_Matches.empty())
  {
      const char *palette[] = {
          "#FF0000", "#FF0100", "#FF0300", "#FF0400", "#FF0500", "#FF0700", "#FF0800", "#FF0900", "#FF0B00", "#FF0C00",
          "#FF0D00", "#FF0F00", "#FF1000", "#FF1100", "#FF1300", "#FF1400", "#FF1500", "#FF1700", "#FF1800", "#FF1900",
          "#FF1B00", "#FF1C00", "#FF1D00", "#FF1F00", "#FF2000", "#FF2100", "#FF2300", "#FF2400", "#FF2500", "#FF2700",
          "#FF2800", "#FF2900", "#FF2B00", "#FF2C00", "#FF2D00", "#FF2F00", "#FF3000", "#FF3100", "#FF3300", "#FF3400",
          "#FF3500", "#FF3700", "#FF3800", "#FF3900", "#FF3B00", "#FF3C00", "#FF3D00", "#FF3F00", "#FF4000", "#FF4100",
          "#FF4300", "#FF4400", "#FF4500", "#FF4700", "#FF4800", "#FF4900", "#FF4B00", "#FF4C00", "#FF4D00", "#FF4F00",
          "#FF5000", "#FF5100", "#FF5300", "#FF5400", "#FF5500", "#FF5700", "#FF5800", "#FF5900", "#FF5B00", "#FF5C00",
          "#FF5D00", "#FF5F00", "#FF6000", "#FF6100", "#FF6300", "#FF6400", "#FF6500", "#FF6700", "#FF6800", "#FF6900",
          "#FF6B00", "#FF6C00", "#FF6D00", "#FF6F00", "#FF7000", "#FF7100", "#FF7300", "#FF7400", "#FF7500", "#FF7700",
          "#FF7800", "#FF7900", "#FF7B00", "#FF7C00", "#FF7D00", "#FF7F00", "#FF8000", "#FF8200", "#FF8300", "#FF8400",
          "#FF8600", "#FF8700", "#FF8800", "#FF8A00", "#FF8B00", "#FF8C00", "#FF8E00", "#FF8F00", "#FF9000", "#FF9200",
          "#FF9300", "#FF9400", "#FF9600", "#FF9700", "#FF9800", "#FF9A00", "#FF9B00", "#FF9C00", "#FF9E00", "#FF9F00",
          "#FFA000", "#FFA200", "#FFA300", "#FFA400", "#FFA600", "#FFA700", "#FFA800", "#FFAA00", "#FFAB00", "#FFAC00",
          "#FFAE00", "#FFAF00", "#FFB000", "#FFB200", "#FFB300", "#FFB400", "#FFB600", "#FFB700", "#FFB800", "#FFBA00",
          "#FFBB00", "#FFBC00", "#FFBE00", "#FFBF00", "#FFC000", "#FFC200", "#FFC300", "#FFC400", "#FFC600", "#FFC700",
          "#FFC800", "#FFCA00", "#FFCB00", "#FFCC00", "#FFCE00", "#FFCF00", "#FFD000", "#FFD200", "#FFD300", "#FFD400",
          "#FFD600", "#FFD700", "#FFD800", "#FFDA00", "#FFDB00", "#FFDC00", "#FFDE00", "#FFDF00", "#FFE000", "#FFE200",
          "#FFE300", "#FFE400", "#FFE600", "#FFE700", "#FFE800", "#FFEA00", "#FFEB00", "#FFEC00", "#FFEE00", "#FFEF00",
          "#FFF000", "#FFF200", "#FFF300", "#FFF400", "#FFF600", "#FFF700", "#FFF800", "#FFFA00", "#FFFB00", "#FFFC00",
          "#FFFE00", "#FFFF00", "#FFFF02", "#FFFF06", "#FFFF0A", "#FFFF0E", "#FFFF12", "#FFFF16", "#FFFF1A", "#FFFF1E",
          "#FFFF22", "#FFFF26", "#FFFF2A", "#FFFF2F", "#FFFF33", "#FFFF37", "#FFFF3B", "#FFFF3F", "#FFFF43", "#FFFF47",
          "#FFFF4B", "#FFFF4F", "#FFFF53", "#FFFF57", "#FFFF5B", "#FFFF5F", "#FFFF63", "#FFFF67", "#FFFF6B", "#FFFF6F",
          "#FFFF73", "#FFFF77", "#FFFF7B", "#FFFF80", "#FFFF84", "#FFFF88", "#FFFF8C", "#FFFF90", "#FFFF94", "#FFFF98",
          "#FFFF9C", "#FFFFA0", "#FFFFA4", "#FFFFA8", "#FFFFAC", "#FFFFB0", "#FFFFB4", "#FFFFB8", "#FFFFBC", "#FFFFC0",
          "#FFFFC4", "#FFFFC8", "#FFFFCC", "#FFFFD0", "#FFFFD5", "#FFFFD9", "#FFFFDD", "#FFFFE1", "#FFFFE5", "#FFFFE9",
          "#FFFFED", "#FFFFF1", "#FFFFF5", "#FFFFF9", "#FFFFFD"
      };
      const size_t n_colors = sizeof(palette) / sizeof(palette[0]);

      size_t min_size = std::numeric_limits<size_t>::max(), max_size = 0;
      for (size_t i = 0; i < NbImages; ++i) {
          for (size_t j = 0; j < NbImages; ++j) {
              matching::PairWiseMatches::const_iterator iter = map_Matches.find(std::make_pair(i, j));
              if (iter != map_Matches.end() && !iter->second.empty()) {
                  const size_t Size = iter->second.size();
                  if (min_size > Size)
                      min_size = Size;
                  if (max_size < Size)
                      max_size = Size;
              }
          }
      }

      float index_scale = static_cast<float>(n_colors) / (max_size - min_size);

    float scaleFactor = 5.0f;
    svgDrawer svgStream((NbImages+3)*5, (NbImages+3)*5);
    // Go along all possible pair
    for (size_t I = 0; I < NbImages; ++I) {
      for (size_t J = 0; J < NbImages; ++J) {
        // If the pair have matches display a blue boxes at I,J position.
        matching::PairWiseMatches::const_iterator iterSearch =
          map_Matches.find(std::make_pair(I,J));
        if (iterSearch != map_Matches.end() && !iterSearch->second.empty())
        {
          // Display as a tooltip: (IndexI, IndexJ NbMatches)
          std::ostringstream os;
          os << "(" << J << "," << I << " " << iterSearch->second.size() <<")";
          size_t color_index = index_scale * (iterSearch->second.size() - min_size) - 1;
          if (color_index >= n_colors) color_index = n_colors - 1;
          svgStream.drawSquare(J*scaleFactor, I*scaleFactor, scaleFactor/2.0f,
            svgStyle().fill(palette[color_index]).noStroke().tooltip(os.str()));
        } // HINT : THINK ABOUT OPACITY [0.4 -> 1.0] TO EXPRESS MATCH COUNT
      }
    }
    // Display axes with 0 -> NbImages annotation : _|
    std::ostringstream osNbImages;   osNbImages << NbImages;
    svgStream.drawText((NbImages+1)*scaleFactor, scaleFactor, scaleFactor, "0", "black");
    svgStream.drawText((NbImages+1)*scaleFactor,
      (NbImages)*scaleFactor - scaleFactor, scaleFactor, osNbImages.str(), "black");
    svgStream.drawLine((NbImages+1)*scaleFactor, 2*scaleFactor,
      (NbImages+1)*scaleFactor, (NbImages)*scaleFactor - 2*scaleFactor,
      svgStyle().stroke("black", 1.0));

    svgStream.drawText(scaleFactor, (NbImages+1)*scaleFactor, scaleFactor, "0", "black");
    svgStream.drawText((NbImages)*scaleFactor - scaleFactor,
      (NbImages+1)*scaleFactor, scaleFactor, osNbImages.str(), "black");
    svgStream.drawLine(2*scaleFactor, (NbImages+1)*scaleFactor,
      (NbImages)*scaleFactor - 2*scaleFactor, (NbImages+1)*scaleFactor,
      svgStyle().stroke("black", 1.0));

    std::ofstream svgFileStream( sOutName.c_str());
    svgFileStream << svgStream.closeSvgFile().str();
  }
}

} // namespace openMVG

#endif // OPENMVG_PAIRWISE_ADJACENCY_DISPLAY_H
