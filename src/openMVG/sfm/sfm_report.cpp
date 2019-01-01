// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm_report.hpp"

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

#include "openMVG/sfm/sfm_data.hpp"

#include "third_party/histogram/histogram.hpp"
#include "third_party/htmlDoc/htmlDoc.hpp"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/vectorGraphics/svgDrawer.hpp"

namespace openMVG {
namespace sfm {

bool Generate_SfM_Report
(
  const SfM_Data & sfm_data,
  const std::string & htmlFilename
)
{
  // Computes:
  // - Statistics per view (mean, max, median residual values)
  // - Global tracks length statistic (length occurence)
  IndexT residualCount = 0;
  Hash_Map<IndexT, std::vector<double>> residuals_per_view;
  std::map<IndexT, IndexT> track_length_occurences;
  for ( const auto & iterTracks : sfm_data.GetLandmarks() )
  {
    const Observations & obs = iterTracks.second.obs;
    track_length_occurences[obs.size()] += 1;
    for ( const auto & itObs : obs )
    {
      const View * view = sfm_data.GetViews().at(itObs.first).get();
      const geometry::Pose3 pose = sfm_data.GetPoseOrDie(view);
      const cameras::IntrinsicBase * intrinsic = sfm_data.GetIntrinsics().at(view->id_intrinsic).get();
      // Use absolute values
      const Vec2 residual = intrinsic->residual(pose(iterTracks.second.X), itObs.second.x).array().abs();
      residuals_per_view[itObs.first].push_back(residual(0));
      residuals_per_view[itObs.first].push_back(residual(1));
      ++residualCount;
    }
  }

  // extract directory from htmlFilename
  const std::string sTableBegin = "<table border=\"1\">",
    sTableEnd = "</table>",
    sRowBegin= "<tr>", sRowEnd = "</tr>",
    sColBegin = "<td>", sColEnd = "</td>",
    sNewLine = "<br>", sFullLine = "<hr>";

  htmlDocument::htmlDocumentStream htmlDocStream("SFM report.");
  htmlDocStream.pushInfo(
  htmlDocument::htmlMarkup("h1", std::string("SFM report.")));
  htmlDocStream.pushInfo(sFullLine);

  htmlDocStream.pushInfo( "Dataset info:" + sNewLine );

  std::ostringstream os;
  os
    << " #views: " << sfm_data.GetViews().size() << sNewLine
    << " #poses: " << sfm_data.GetPoses().size() << sNewLine
    << " #intrinsics: " << sfm_data.GetIntrinsics().size() << sNewLine
    << " #tracks: " << sfm_data.GetLandmarks().size() << sNewLine
    << " #residuals: " << residualCount << sNewLine;

  htmlDocStream.pushInfo( os.str() );
  htmlDocStream.pushInfo( sFullLine );

  // Track length statistics
  if (!track_length_occurences.empty() && !sfm_data.GetLandmarks().empty())
  {
    const IndexT min_track_length = track_length_occurences.cbegin()->first;
    const IndexT max_track_length = track_length_occurences.crbegin()->first;
    // Compute the mean track length
    IndexT sum = 0;
    for (const auto & pair_length_occurence : track_length_occurences)
    {
      sum += pair_length_occurence.first * pair_length_occurence.second;
    }
    const IndexT mean_track_length = sum / static_cast<double>(sfm_data.GetLandmarks().size());

    htmlDocStream.pushInfo( "Track length statistics:" );
    htmlDocStream.pushInfo(sTableBegin);
    os.str("");
    os
      << sRowBegin
      << sColBegin + "min" + sColEnd
      << sColBegin + "mean" + sColEnd
      << sColBegin + "max" + sColEnd
      << sRowEnd
      << sRowBegin
      << sColBegin << min_track_length << sColEnd
      << sColBegin << mean_track_length << sColEnd
      << sColBegin << max_track_length << sColEnd
      << sRowEnd
      << sTableEnd;
    htmlDocStream.pushInfo( os.str() );
  }

  htmlDocStream.pushInfo( sFullLine );

  // Residuals statistics and observation count per view
  htmlDocStream.pushInfo( sTableBegin);
  os.str("");
  os << sRowBegin
    << sColBegin + "IdView" + sColEnd
    << sColBegin + "Basename" + sColEnd
    << sColBegin + "#Observations" + sColEnd
    << sColBegin + "Residuals min" + sColEnd
    << sColBegin + "Residuals median" + sColEnd
    << sColBegin + "Residuals mean" + sColEnd
    << sColBegin + "Residuals max" + sColEnd
    << sRowEnd;
  htmlDocStream.pushInfo( os.str() );

  for (const auto & iterV : sfm_data.GetViews())
  {
    const View * v = iterV.second.get();
    const IndexT id_view = v->id_view;

    os.str("");
    os << sRowBegin
      << sColBegin << id_view << sColEnd
      << sColBegin + stlplus::basename_part(v->s_Img_path) + sColEnd;

    // IdView | basename | #Observations | residuals min | residual median | residual max
    if (sfm_data.IsPoseAndIntrinsicDefined(v))
    {
      if (residuals_per_view.find(id_view) != residuals_per_view.end())
      {
        const std::vector<double> & residuals = residuals_per_view.at(id_view);
        if (!residuals.empty())
        {
          double min, max, mean, median;
          minMaxMeanMedian(residuals.begin(), residuals.end(), min, max, mean, median);
          os << sColBegin << residuals.size()/2 << sColEnd // #observations
            << sColBegin << min << sColEnd
            << sColBegin << median << sColEnd
            << sColBegin << mean << sColEnd
            << sColBegin << max <<sColEnd;
        }
      }
    }
    os << sRowEnd;
    htmlDocStream.pushInfo( os.str() );
  }
  htmlDocStream.pushInfo( sTableEnd );
  htmlDocStream.pushInfo( sFullLine );

  // combine all residual values into one vector
  // export the SVG histogram
  {
    IndexT residualCount = 0;
    for (const auto & residual_per_view_it : residuals_per_view)
    {
      residualCount += residual_per_view_it.second.size();
    }
    // Concat per view residual values into one vector
    std::vector<double> residuals(residualCount);
    residualCount = 0;
    for (const auto & residual_per_view_it : residuals_per_view)
    {
      std::copy(residual_per_view_it.second.begin(),
        residual_per_view_it.second.end(),
        residuals.begin() + residualCount);
      residualCount += residual_per_view_it.second.size();
    }
    if (!residuals.empty())
    {
      // RMSE computation
      const Eigen::Map<Eigen::RowVectorXd> residuals_mapping(&residuals[0], residuals.size());
      const double RMSE = std::sqrt(residuals_mapping.squaredNorm() / (double)residuals.size());
      os.str("");
      os << sFullLine << "SfM Scene RMSE: " << RMSE << sFullLine;
      htmlDocStream.pushInfo(os.str());

      const double maxRange = *max_element(residuals.cbegin(), residuals.cend());
      Histogram<double> histo(0.0, maxRange, 100);
      histo.Add(residuals.cbegin(), residuals.cend());

      svg::svgHisto svg_Histo;
      svg_Histo.draw(histo.GetHist(), std::pair<float,float>(0.f, maxRange),
        stlplus::create_filespec(stlplus::folder_part(htmlFilename), "residuals_histogram", "svg"),
        600, 200);

      os.str("");
      os << sNewLine<< "Residuals histogram" << sNewLine;
      os << "<img src=\""
        << "residuals_histogram.svg"
        << "\" height=\"300\" width =\"800\">\n";
      htmlDocStream.pushInfo(os.str());
    }
  }

  std::ofstream htmlFileStream(htmlFilename);
  htmlFileStream << htmlDocStream.getDoc();
  return !htmlFileStream.good();
}

} // namespace sfm
} // namespace openMVG
