// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm_report.hpp"

#include "openMVG/numeric/numeric.h"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/system/logger.hpp"
#include "openMVG/vector_graphics/svgHistogram.hpp"

#include "third_party/histogram/histogram.hpp"
#include "third_party/htmlDoc/htmlDoc.hpp"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace openMVG {
namespace sfm {

namespace {

/**
 * @brief Compute linear-interpolated percentile (0 to 100) from a sorted vector
 */
template <typename T>
double ComputePercentile(const std::vector<T>& sorted_vec, double pct)
{
  if (sorted_vec.empty()) return 0.0;
  if (sorted_vec.size() == 1) return static_cast<double>(sorted_vec[0]);
  const double rank = (pct / 100.0) * static_cast<double>(sorted_vec.size() - 1);
  const size_t lower = static_cast<size_t>(std::floor(rank));
  const size_t upper = static_cast<size_t>(std::ceil(rank));
  const double weight = rank - static_cast<double>(lower);
  return (1.0 - weight) * static_cast<double>(sorted_vec[lower]) + weight * static_cast<double>(sorted_vec[upper]);
}

} // namespace

bool Generate_SfM_Report
(
  const SfM_Data & sfm_data,
  const std::string & htmlFilename
)
{
  // Computes:
  // - Statistics per view (mean, max, median residual values)
  // - Detailed global track length statistics (mean, median, percentiles, distribution, SVG)
  // - Parallax angles for triangulation quality evaluation
  // - Camera covisibility graph density
  // - 1D and 2D reprojection residuals (RMSE, mean, median, percentiles, error bands)

  if (sfm_data.GetViews().empty())
  {
    OPENMVG_LOG_WARNING << "Cannot generate SfM report: scene has no views.";
    return false;
  }

  IndexT residualCount = 0;
  Hash_Map<IndexT, std::vector<double>> residuals_per_view;
  std::vector<double> residuals_2d_all;
  std::map<IndexT, IndexT> track_length_occurrences;
  std::vector<double> all_track_lengths;
  all_track_lengths.reserve(sfm_data.GetLandmarks().size());
  IndexT total_observations = 0;

  // Covisibility tracking: map of camera pair (c1, c2) with c1 < c2 -> shared landmarks count
  std::map<std::pair<IndexT, IndexT>, IndexT> covisibility_counts;

  // Parallax angles for triangulation
  std::vector<double> max_parallax_angles;
  max_parallax_angles.reserve(sfm_data.GetLandmarks().size());

  for ( const auto & iterTracks : sfm_data.GetLandmarks() )
  {
    const Observations & obs = iterTracks.second.obs;
    const size_t track_len = obs.size();
    track_length_occurrences[track_len] += 1;
    all_track_lengths.push_back(static_cast<double>(track_len));
    total_observations += track_len;

    // Track covisibility between camera pairs
    std::vector<IndexT> obs_view_ids;
    obs_view_ids.reserve(track_len);

    // Collect ray directions for parallax angle calculation
    std::vector<Vec3> ray_dirs;
    ray_dirs.reserve(track_len);

    for ( const auto & itObs : obs )
    {
      const auto itView = sfm_data.GetViews().find(itObs.first);
      if (itView == sfm_data.GetViews().end())
        continue;
      const View * view = itView->second.get();
      if (!sfm_data.IsPoseAndIntrinsicDefined(view))
        continue;

      obs_view_ids.push_back(itObs.first);

      const geometry::Pose3 pose = sfm_data.GetPoseOrDie(view);
      const cameras::IntrinsicBase * intrinsic = sfm_data.GetIntrinsics().at(view->id_intrinsic).get();

      // Residual vector
      const Vec2 residual = intrinsic->residual(pose(iterTracks.second.X), itObs.second.x);
      const Vec2 abs_residual = residual.array().abs();
      residuals_per_view[itObs.first].push_back(abs_residual(0));
      residuals_per_view[itObs.first].push_back(abs_residual(1));
      residuals_2d_all.push_back(residual.norm());
      ++residualCount;

      // Optical ray for parallax calculation
      const Vec3 ray = (iterTracks.second.X - pose.center());
      const double norm = ray.norm();
      if (norm > 1e-8)
      {
        ray_dirs.push_back(ray / norm);
      }
    }

    // Accumulate covisibility for pairs
    for (size_t i = 0; i < obs_view_ids.size(); ++i)
    {
      for (size_t j = i + 1; j < obs_view_ids.size(); ++j)
      {
        const IndexT id1 = std::min(obs_view_ids[i], obs_view_ids[j]);
        const IndexT id2 = std::max(obs_view_ids[i], obs_view_ids[j]);
        covisibility_counts[std::make_pair(id1, id2)] += 1;
      }
    }

    // Parallax angle for this track
    if (ray_dirs.size() >= 2)
    {
      double max_ang_deg = 0.0;
      for (size_t i = 0; i < ray_dirs.size(); ++i)
      {
        for (size_t j = i + 1; j < ray_dirs.size(); ++j)
        {
          const double dot = std::max(-1.0, std::min(1.0, ray_dirs[i].dot(ray_dirs[j])));
          const double ang = std::acos(dot) * (180.0 / M_PI);
          if (ang > max_ang_deg)
            max_ang_deg = ang;
        }
      }
      max_parallax_angles.push_back(max_ang_deg);
    }
  }

  // HTML layout tags
  const std::string sTableBegin = "<table border=\"1\" cellpadding=\"4\" cellspacing=\"0\" style=\"border-collapse:collapse; margin-bottom:15px;\">",
    sTableEnd = "</table>",
    sRowBegin= "<tr>", sRowEnd = "</tr>",
    sColBegin = "<td>", sColEnd = "</td>",
    sColHeadBegin = "<th style=\"background-color:#f2f2f2;\">", sColHeadEnd = "</th>",
    sNewLine = "<br>", sFullLine = "<hr>";

  htmlDocument::htmlDocumentStream htmlDocStream("SFM report.");
  htmlDocStream.pushInfo(htmlDocument::htmlMarkup("h1", std::string("SFM report.")));
  htmlDocStream.pushInfo(sFullLine);

  htmlDocStream.pushInfo("<h2>Dataset info:</h2>" + sNewLine);

  std::ostringstream os;
  os
    << " #views: " << sfm_data.GetViews().size() << sNewLine
    << " #poses: " << sfm_data.GetPoses().size() << sNewLine
    << " #intrinsics: " << sfm_data.GetIntrinsics().size() << sNewLine
    << " #tracks (3D points): " << sfm_data.GetLandmarks().size() << sNewLine
    << " #observations (2D features): " << total_observations << sNewLine
    << " #residuals: " << residualCount << sNewLine;

  htmlDocStream.pushInfo( os.str() );
  htmlDocStream.pushInfo( sFullLine );

  // --------------------------------------------------------------------------
  // 1. Detailed Track Length Statistics
  // --------------------------------------------------------------------------
  double mean_track_length = 0.0;
  double median_track_length = 0.0;
  double min_track_length = 0.0;
  double max_track_length = 0.0;
  double p25_track = 0.0, p50_track = 0.0, p75_track = 0.0, p90_track = 0.0, p95_track = 0.0, p99_track = 0.0;
  const size_t num_tracks = all_track_lengths.size();

  struct CategoryInfo {
    std::string label;
    IndexT count;
    IndexT obs;
    CategoryInfo(const std::string& l = "", IndexT c = 0, IndexT o = 0)
      : label(l), count(c), obs(o) {}
  };
  CategoryInfo cat2   { "2 views (minimal)", 0, 0 };
  CategoryInfo cat3   { "3 views (short)", 0, 0 };
  CategoryInfo cat4_5 { "4-5 views (medium)", 0, 0 };
  CategoryInfo cat6_9 { "6-9 views (good multi-view)", 0, 0 };
  CategoryInfo cat10p { ">= 10 views (robust)", 0, 0 };

  if (!all_track_lengths.empty() && num_tracks > 0)
  {
    std::sort(all_track_lengths.begin(), all_track_lengths.end());
    min_track_length = all_track_lengths.front();
    max_track_length = all_track_lengths.back();
    mean_track_length = static_cast<double>(total_observations) / static_cast<double>(num_tracks);
    median_track_length = ComputePercentile(all_track_lengths, 50.0);

    double sum_sq_diff = 0.0;
    for (const double len : all_track_lengths)
    {
      const double diff = len - mean_track_length;
      sum_sq_diff += diff * diff;
    }
    const double stddev_track_length = std::sqrt(sum_sq_diff / static_cast<double>(num_tracks));

    p25_track = ComputePercentile(all_track_lengths, 25.0);
    p50_track = ComputePercentile(all_track_lengths, 50.0);
    p75_track = ComputePercentile(all_track_lengths, 75.0);
    p90_track = ComputePercentile(all_track_lengths, 90.0);
    p95_track = ComputePercentile(all_track_lengths, 95.0);
    p99_track = ComputePercentile(all_track_lengths, 99.0);

    htmlDocStream.pushInfo( "<h2>Track length statistics:</h2>" );

    // Summary table
    htmlDocStream.pushInfo( "<h3>Summary</h3>" );
    htmlDocStream.pushInfo(sTableBegin);
    os.str("");
    os << std::fixed << std::setprecision(2);
    os
      << sRowBegin
      << sColHeadBegin + "min" + sColHeadEnd
      << sColHeadBegin + "mean" + sColHeadEnd
      << sColHeadBegin + "median" + sColHeadEnd
      << sColHeadBegin + "stddev" + sColHeadEnd
      << sColHeadBegin + "P25" + sColHeadEnd
      << sColHeadBegin + "P50" + sColHeadEnd
      << sColHeadBegin + "P75" + sColHeadEnd
      << sColHeadBegin + "P90" + sColHeadEnd
      << sColHeadBegin + "P95" + sColHeadEnd
      << sColHeadBegin + "P99" + sColHeadEnd
      << sColHeadBegin + "max" + sColHeadEnd
      << sRowEnd
      << sRowBegin
      << sColBegin << static_cast<IndexT>(min_track_length) << sColEnd
      << sColBegin << mean_track_length << sColEnd
      << sColBegin << median_track_length << sColEnd
      << sColBegin << stddev_track_length << sColEnd
      << sColBegin << p25_track << sColEnd
      << sColBegin << p50_track << sColEnd
      << sColBegin << p75_track << sColEnd
      << sColBegin << p90_track << sColEnd
      << sColBegin << p95_track << sColEnd
      << sColBegin << p99_track << sColEnd
      << sColBegin << static_cast<IndexT>(max_track_length) << sColEnd
      << sRowEnd
      << sTableEnd;
    htmlDocStream.pushInfo( os.str() );

    // Categorized breakdown table
    for (const auto & pair : track_length_occurrences)
    {
      const IndexT len = pair.first;
      const IndexT cnt = pair.second;
      const IndexT o   = len * cnt;
      if (len == 2) { cat2.count += cnt; cat2.obs += o; }
      else if (len == 3) { cat3.count += cnt; cat3.obs += o; }
      else if (len <= 5) { cat4_5.count += cnt; cat4_5.obs += o; }
      else if (len <= 9) { cat6_9.count += cnt; cat6_9.obs += o; }
      else { cat10p.count += cnt; cat10p.obs += o; }
    }

    const CategoryInfo all_cats[5] = { cat2, cat3, cat4_5, cat6_9, cat10p };

    htmlDocStream.pushInfo( "<h3>Track Length Categories</h3>" );
    htmlDocStream.pushInfo(sTableBegin);
    os.str("");
    os
      << sRowBegin
      << sColHeadBegin + "Category" + sColHeadEnd
      << sColHeadBegin + "#Tracks" + sColHeadEnd
      << sColHeadBegin + "% of Tracks" + sColHeadEnd
      << sColHeadBegin + "#Observations" + sColHeadEnd
      << sColHeadBegin + "% of Observations" + sColHeadEnd
      << sRowEnd;
    for (const auto & cat : all_cats)
    {
      const double pct_tracks = 100.0 * static_cast<double>(cat.count) / static_cast<double>(num_tracks);
      const double pct_obs = (total_observations > 0)
        ? 100.0 * static_cast<double>(cat.obs) / static_cast<double>(total_observations) : 0.0;
      os
        << sRowBegin
        << sColBegin + cat.label + sColEnd
        << sColBegin << cat.count << sColEnd
        << sColBegin << pct_tracks << " %" + sColEnd
        << sColBegin << cat.obs << sColEnd
        << sColBegin << pct_obs << " %" + sColEnd
        << sRowEnd;
    }
    os << sTableEnd;
    htmlDocStream.pushInfo( os.str() );

    // Detailed occurrence table per track length
    htmlDocStream.pushInfo( "<h3>Detailed Distribution per Length</h3>" );
    htmlDocStream.pushInfo(sTableBegin);
    os.str("");
    os
      << sRowBegin
      << sColHeadBegin + "Length" + sColHeadEnd
      << sColHeadBegin + "#Tracks" + sColHeadEnd
      << sColHeadBegin + "% of Tracks" + sColHeadEnd
      << sColHeadBegin + "Cumulative %" + sColHeadEnd
      << sColHeadBegin + "#Observations" + sColHeadEnd
      << sColHeadBegin + "% of Observations" + sColHeadEnd
      << sRowEnd;

    IndexT cumul_tracks = 0;
    for (const auto & pair : track_length_occurrences)
    {
      const IndexT len = pair.first;
      const IndexT cnt = pair.second;
      cumul_tracks += cnt;
      const IndexT o = len * cnt;
      const double pct_tracks = 100.0 * static_cast<double>(cnt) / static_cast<double>(num_tracks);
      const double cumul_pct  = 100.0 * static_cast<double>(cumul_tracks) / static_cast<double>(num_tracks);
      const double pct_obs    = (total_observations > 0)
        ? 100.0 * static_cast<double>(o) / static_cast<double>(total_observations) : 0.0;

      os
        << sRowBegin
        << sColBegin << len << sColEnd
        << sColBegin << cnt << sColEnd
        << sColBegin << pct_tracks << " %" + sColEnd
        << sColBegin << cumul_pct << " %" + sColEnd
        << sColBegin << o << sColEnd
        << sColBegin << pct_obs << " %" + sColEnd
        << sRowEnd;
    }
    os << sTableEnd;
    htmlDocStream.pushInfo( os.str() );

    // Generate tracks histogram SVG
    if (max_track_length >= min_track_length)
    {
      const size_t num_bins = static_cast<size_t>(max_track_length - min_track_length + 1);
      std::vector<IndexT> track_hist(num_bins, 0);
      for (const auto & pair : track_length_occurrences)
      {
        track_hist[pair.first - static_cast<IndexT>(min_track_length)] = pair.second;
      }
      svg::drawHistogram(track_hist,
        std::make_pair(static_cast<float>(min_track_length), static_cast<float>(max_track_length)),
        stlplus::create_filespec(stlplus::folder_part(htmlFilename), "tracks_histogram", "svg"),
        600, 200);

      os.str("");
      os << sNewLine << "<b>Tracks length histogram:</b>" << sNewLine;
      os << "<img src=\"tracks_histogram.svg\" height=\"250\" width=\"700\">\n";
      htmlDocStream.pushInfo(os.str());
    }
  }

  htmlDocStream.pushInfo( sFullLine );

  // --------------------------------------------------------------------------
  // 2. Triangulation Parallax Angles Statistics
  // --------------------------------------------------------------------------
  double min_ang = 0.0, max_ang = 0.0, mean_ang = 0.0, median_ang = 0.0;
  const size_t N_ang = max_parallax_angles.size();
  IndexT ang_lt2 = 0, ang_2_5 = 0, ang_5_10 = 0, ang_10_20 = 0, ang_ge20 = 0;

  if (!max_parallax_angles.empty())
  {
    std::sort(max_parallax_angles.begin(), max_parallax_angles.end());
    minMaxMeanMedian(max_parallax_angles.begin(), max_parallax_angles.end(), min_ang, max_ang, mean_ang, median_ang);

    double sum_sq_diff_ang = 0.0;
    for (const double ang : max_parallax_angles)
    {
      const double diff = ang - mean_ang;
      sum_sq_diff_ang += diff * diff;
    }
    const double stddev_ang = std::sqrt(sum_sq_diff_ang / static_cast<double>(N_ang));

    const double p25_ang = ComputePercentile(max_parallax_angles, 25.0);
    const double p50_ang = ComputePercentile(max_parallax_angles, 50.0);
    const double p75_ang = ComputePercentile(max_parallax_angles, 75.0);
    const double p90_ang = ComputePercentile(max_parallax_angles, 90.0);
    const double p95_ang = ComputePercentile(max_parallax_angles, 95.0);

    htmlDocStream.pushInfo( "<h2>Triangulation Parallax Angles (&theta;<sub>max</sub>):</h2>" );
    htmlDocStream.pushInfo( "<h3>Summary (degrees)</h3>" );
    htmlDocStream.pushInfo(sTableBegin);
    os.str("");
    os << std::fixed << std::setprecision(2);
    os
      << sRowBegin
      << sColHeadBegin + "min" + sColHeadEnd
      << sColHeadBegin + "mean" + sColHeadEnd
      << sColHeadBegin + "median" + sColHeadEnd
      << sColHeadBegin + "stddev" + sColHeadEnd
      << sColHeadBegin + "P25" + sColHeadEnd
      << sColHeadBegin + "P50" + sColHeadEnd
      << sColHeadBegin + "P75" + sColHeadEnd
      << sColHeadBegin + "P90" + sColHeadEnd
      << sColHeadBegin + "P95" + sColHeadEnd
      << sColHeadBegin + "max" + sColHeadEnd
      << sRowEnd
      << sRowBegin
      << sColBegin << min_ang << "&deg;" + sColEnd
      << sColBegin << mean_ang << "&deg;" + sColEnd
      << sColBegin << median_ang << "&deg;" + sColEnd
      << sColBegin << stddev_ang << "&deg;" + sColEnd
      << sColBegin << p25_ang << "&deg;" + sColEnd
      << sColBegin << p50_ang << "&deg;" + sColEnd
      << sColBegin << p75_ang << "&deg;" + sColEnd
      << sColBegin << p90_ang << "&deg;" + sColEnd
      << sColBegin << p95_ang << "&deg;" + sColEnd
      << sColBegin << max_ang << "&deg;" + sColEnd
      << sRowEnd
      << sTableEnd;
    htmlDocStream.pushInfo( os.str() );

    // Parallax angle bands
    for (const double a : max_parallax_angles)
    {
      if (a < 2.0) ++ang_lt2;
      else if (a < 5.0) ++ang_2_5;
      else if (a < 10.0) ++ang_5_10;
      else if (a < 20.0) ++ang_10_20;
      else ++ang_ge20;
    }

    htmlDocStream.pushInfo( "<h3>Parallax Angle Distribution</h3>" );
    htmlDocStream.pushInfo(sTableBegin);
    os.str("");
    auto add_ang_row = [&](const std::string& band, IndexT cnt) {
      const double pct = 100.0 * static_cast<double>(cnt) / static_cast<double>(N_ang);
      os << sRowBegin
         << sColBegin + band + sColEnd
         << sColBegin << cnt << sColEnd
         << sColBegin << pct << " %" + sColEnd
         << sRowEnd;
    };
    os
      << sRowBegin
      << sColHeadBegin + "Angle Band" + sColHeadEnd
      << sColHeadBegin + "#Tracks" + sColHeadEnd
      << sColHeadBegin + "% of Tracks" + sColHeadEnd
      << sRowEnd;
    add_ang_row("< 2.0&deg; (fragile / small baseline)", ang_lt2);
    add_ang_row("[2.0&deg;, 5.0&deg;)", ang_2_5);
    add_ang_row("[5.0&deg;, 10.0&deg;)", ang_5_10);
    add_ang_row("[10.0&deg;, 20.0&deg;) (well conditioned)", ang_10_20);
    add_ang_row(">= 20.0&deg; (wide baseline)", ang_ge20);
    os << sTableEnd;
    htmlDocStream.pushInfo( os.str() );

    // Parallax angles histogram SVG
    if (max_ang > min_ang)
    {
      Histogram<double> histo_ang(0.0, max_ang, 50);
      histo_ang.Add(max_parallax_angles.cbegin(), max_parallax_angles.cend());
      svg::drawHistogram(histo_ang.GetHist(), std::pair<float,float>(0.f, static_cast<float>(max_ang)),
        stlplus::create_filespec(stlplus::folder_part(htmlFilename), "parallax_angles_histogram", "svg"),
        600, 200);

      os.str("");
      os << sNewLine << "<b>Parallax angles histogram (degrees):</b>" << sNewLine;
      os << "<img src=\"parallax_angles_histogram.svg\" height=\"250\" width=\"700\">\n";
      htmlDocStream.pushInfo(os.str());
    }

    htmlDocStream.pushInfo( sFullLine );
  }

  // --------------------------------------------------------------------------
  // 3. Camera Covisibility Graph Statistics
  // --------------------------------------------------------------------------
  if (!covisibility_counts.empty())
  {
    std::vector<IndexT> shared_pts;
    shared_pts.reserve(covisibility_counts.size());
    IndexT covis_ge30 = 0;
    for (const auto & pair : covisibility_counts)
    {
      shared_pts.push_back(pair.second);
      if (pair.second >= 30)
        ++covis_ge30;
    }
    std::sort(shared_pts.begin(), shared_pts.end());
    const size_t N_pairs = shared_pts.size();
    const double mean_shared = static_cast<double>(std::accumulate(shared_pts.begin(), shared_pts.end(), IndexT(0))) / static_cast<double>(N_pairs);
    const double median_shared = ComputePercentile(shared_pts, 50.0);

    htmlDocStream.pushInfo( "<h2>Camera Covisibility Graph:</h2>" );
    htmlDocStream.pushInfo(sTableBegin);
    os.str("");
    os << std::fixed << std::setprecision(2);
    os
      << sRowBegin
      << sColHeadBegin + "Metric" + sColHeadEnd
      << sColHeadBegin + "Value" + sColHeadEnd
      << sRowEnd
      << sRowBegin
      << sColBegin + "Covisible camera pairs (sharing >= 1 point)" + sColEnd
      << sColBegin << N_pairs << sColEnd
      << sRowEnd
      << sRowBegin
      << sColBegin + "Covisible camera pairs sharing >= 30 points" + sColEnd
      << sColBegin << covis_ge30 << sColEnd
      << sRowEnd
      << sRowBegin
      << sColBegin + "Shared points per covisible pair (min / median / mean / max)" + sColEnd
      << sColBegin << shared_pts.front() << " / " << median_shared << " / " << mean_shared << " / " << shared_pts.back() << sColEnd
      << sRowEnd
      << sTableEnd;
    htmlDocStream.pushInfo( os.str() );
    htmlDocStream.pushInfo( sFullLine );
  }

  // --------------------------------------------------------------------------
  // 4. Residuals Statistics per View
  // --------------------------------------------------------------------------
  htmlDocStream.pushInfo("<h2>Residuals Statistics per View:</h2>");
  htmlDocStream.pushInfo( sTableBegin);
  os.str("");
  os << sRowBegin
    << sColHeadBegin + "IdView" + sColHeadEnd
    << sColHeadBegin + "Basename" + sColHeadEnd
    << sColHeadBegin + "#Observations" + sColHeadEnd
    << sColHeadBegin + "Residuals min" + sColHeadEnd
    << sColHeadBegin + "Residuals median" + sColHeadEnd
    << sColHeadBegin + "Residuals mean" + sColHeadEnd
    << sColHeadBegin + "Residuals max" + sColHeadEnd
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

  // --------------------------------------------------------------------------
  // 5. Global Residuals & 2D Reprojection Error Analysis
  // --------------------------------------------------------------------------
  double rmse_2d = 0.0, mean_e = 0.0, median_e = 0.0, max_e = 0.0;
  double p90_e = 0.0, p95_e = 0.0, p99_e = 0.0;
  IndexT b_gt2 = 0;
  const size_t N_res = residuals_2d_all.size();

  if (!residuals_2d_all.empty())
  {
    std::sort(residuals_2d_all.begin(), residuals_2d_all.end());
    double sum_sq_2d = 0.0;
    for (const double e : residuals_2d_all)
      sum_sq_2d += e * e;
    rmse_2d = std::sqrt(sum_sq_2d / static_cast<double>(N_res));

    double min_e;
    minMaxMeanMedian(residuals_2d_all.begin(), residuals_2d_all.end(), min_e, max_e, mean_e, median_e);

    const double p25_e = ComputePercentile(residuals_2d_all, 25.0);
    const double p50_e = ComputePercentile(residuals_2d_all, 50.0);
    const double p75_e = ComputePercentile(residuals_2d_all, 75.0);
    p90_e = ComputePercentile(residuals_2d_all, 90.0);
    p95_e = ComputePercentile(residuals_2d_all, 95.0);
    p99_e = ComputePercentile(residuals_2d_all, 99.0);

    IndexT b05 = 0, b10 = 0, b20 = 0;
    for (const double e : residuals_2d_all)
    {
      if (e <= 0.5) ++b05;
      else if (e <= 1.0) ++b10;
      else if (e <= 2.0) ++b20;
      else ++b_gt2;
    }

    htmlDocStream.pushInfo( "<h2>Reprojection Error (2D Norm) Statistics:</h2>" );
    htmlDocStream.pushInfo( "<h3>Summary (pixels)</h3>" );
    htmlDocStream.pushInfo(sTableBegin);
    os.str("");
    os << std::fixed << std::setprecision(4);
    os
      << sRowBegin
      << sColHeadBegin + "RMSE (2D)" + sColHeadEnd
      << sColHeadBegin + "mean" + sColHeadEnd
      << sColHeadBegin + "median" + sColHeadEnd
      << sColHeadBegin + "P25" + sColHeadEnd
      << sColHeadBegin + "P50" + sColHeadEnd
      << sColHeadBegin + "P75" + sColHeadEnd
      << sColHeadBegin + "P90" + sColHeadEnd
      << sColHeadBegin + "P95" + sColHeadEnd
      << sColHeadBegin + "P99" + sColHeadEnd
      << sColHeadBegin + "max" + sColHeadEnd
      << sRowEnd
      << sRowBegin
      << sColBegin << rmse_2d << " px" + sColEnd
      << sColBegin << mean_e << " px" + sColEnd
      << sColBegin << median_e << " px" + sColEnd
      << sColBegin << p25_e << " px" + sColEnd
      << sColBegin << p50_e << " px" + sColEnd
      << sColBegin << p75_e << " px" + sColEnd
      << sColBegin << p90_e << " px" + sColEnd
      << sColBegin << p95_e << " px" + sColEnd
      << sColBegin << p99_e << " px" + sColEnd
      << sColBegin << max_e << " px" + sColEnd
      << sRowEnd
      << sTableEnd;
    htmlDocStream.pushInfo( os.str() );

    // Error bands table
    htmlDocStream.pushInfo( "<h3>Reprojection Error Distribution</h3>" );
    htmlDocStream.pushInfo(sTableBegin);
    os.str("");
    os << std::fixed << std::setprecision(2);
    auto add_err_row = [&](const std::string& band, IndexT cnt) {
      const double pct = 100.0 * static_cast<double>(cnt) / static_cast<double>(N_res);
      os << sRowBegin
         << sColBegin + band + sColEnd
         << sColBegin << cnt << sColEnd
         << sColBegin << pct << " %" + sColEnd
         << sRowEnd;
    };
    os
      << sRowBegin
      << sColHeadBegin + "Error Range" + sColHeadEnd
      << sColHeadBegin + "#Observations" + sColHeadEnd
      << sColHeadBegin + "% of Observations" + sColHeadEnd
      << sRowEnd;
    add_err_row("[0.0, 0.5] px (subpixel fine)", b05);
    add_err_row("(0.5, 1.0] px (good)", b10);
    add_err_row("(1.0, 2.0] px (acceptable)", b20);
    add_err_row("> 2.0 px (high residual)", b_gt2);
    os << sTableEnd;
    htmlDocStream.pushInfo( os.str() );
  }

  // Combine 1D residual values for historical compatibility and 1D RMSE SVG histogram
  {
    IndexT residual1DCount = 0;
    for (const auto & residual_per_view_it : residuals_per_view)
    {
      residual1DCount += residual_per_view_it.second.size();
    }
    std::vector<double> residuals(residual1DCount);
    residual1DCount = 0;
    for (const auto & residual_per_view_it : residuals_per_view)
    {
      std::copy(residual_per_view_it.second.begin(),
        residual_per_view_it.second.end(),
        residuals.begin() + residual1DCount);
      residual1DCount += residual_per_view_it.second.size();
    }
    if (!residuals.empty())
    {
      // Historical 1D RMSE computation
      const Eigen::Map<Eigen::RowVectorXd> residuals_mapping(&residuals[0], residuals.size());
      const double RMSE = std::sqrt(residuals_mapping.squaredNorm() / (double)residuals.size());
      os.str("");
      os << sFullLine << "<b>SfM Scene RMSE (1D coordinate residual):</b> " << std::fixed << std::setprecision(4) << RMSE << sNewLine;
      htmlDocStream.pushInfo(os.str());

      const double maxRange = *max_element(residuals.cbegin(), residuals.cend());
      Histogram<double> histo(0.0, maxRange, 100);
      histo.Add(residuals.cbegin(), residuals.cend());

      svg::drawHistogram(histo.GetHist(), std::pair<float,float>(0.f, static_cast<float>(maxRange)),
        stlplus::create_filespec(stlplus::folder_part(htmlFilename), "residuals_histogram", "svg"),
        600, 200);

      os.str("");
      os << sNewLine << "<b>Residuals histogram:</b>" << sNewLine;
      os << "<img src=\""
        << "residuals_histogram.svg"
        << "\" height=\"300\" width =\"800\">\n";
      htmlDocStream.pushInfo(os.str());
    }
  }

  std::ofstream htmlFileStream(htmlFilename);
  if (!htmlFileStream.is_open())
  {
    OPENMVG_LOG_ERROR << "Failed to open HTML report file for writing: " << htmlFilename;
    return false;
  }
  htmlFileStream << htmlDocStream.getDoc();

  // --------------------------------------------------------------------------
  // Console logging (terminal summary)
  // --------------------------------------------------------------------------
  OPENMVG_LOG_INFO
    << "=================================================================\n"
    << "SfM Reconstruction Report Summary:\n"
    << "-----------------------------------------------------------------\n"
    << " - Views: " << sfm_data.GetViews().size() << ", Poses: " << sfm_data.GetPoses().size() << "\n"
    << " - Landmarks (3D points): " << sfm_data.GetLandmarks().size() << "\n"
    << " - Observations (2D features): " << total_observations << "\n"
    << " - Track length: min " << static_cast<IndexT>(min_track_length)
    << ", mean " << std::fixed << std::setprecision(2) << mean_track_length
    << ", median " << median_track_length
    << ", max " << static_cast<IndexT>(max_track_length) << "\n"
    << "   P25: " << p25_track << ", P50: " << p50_track << ", P75: " << p75_track
    << ", P90: " << p90_track << ", P95: " << p95_track << ", P99: " << p99_track << "\n"
    << "   Categories: 2-views: " << cat2.count << " (" << (num_tracks > 0 ? 100.0 * cat2.count / num_tracks : 0.0) << "%), "
    << "3-views: " << cat3.count << " (" << (num_tracks > 0 ? 100.0 * cat3.count / num_tracks : 0.0) << "%), "
    << "4-5 views: " << cat4_5.count << " (" << (num_tracks > 0 ? 100.0 * cat4_5.count / num_tracks : 0.0) << "%), "
    << "6-9 views: " << cat6_9.count << " (" << (num_tracks > 0 ? 100.0 * cat6_9.count / num_tracks : 0.0) << "%), "
    << "10+ views: " << cat10p.count << " (" << (num_tracks > 0 ? 100.0 * cat10p.count / num_tracks : 0.0) << "%)\n";

  if (!residuals_2d_all.empty())
  {
    OPENMVG_LOG_INFO
      << " - 2D Reprojection RMSE: " << std::fixed << std::setprecision(4) << rmse_2d << " px\n"
      << "   mean: " << mean_e << " px, median: " << median_e << " px\n"
      << "   P90: " << p90_e << " px, P95: " << p95_e << " px, P99: " << p99_e << " px, max: " << max_e << " px\n"
      << "   Residuals > 2.0 px: " << b_gt2 << " (" << std::setprecision(2) << 100.0 * b_gt2 / N_res << "%)\n";
  }
  if (!max_parallax_angles.empty())
  {
    OPENMVG_LOG_INFO
      << " - Parallax angles: mean " << std::fixed << std::setprecision(2) << mean_ang << " deg"
      << ", median " << median_ang << " deg"
      << ", max " << max_ang << " deg"
      << " (>= 10 deg: " << (N_ang > 0 ? 100.0 * (ang_10_20 + ang_ge20) / N_ang : 0.0) << "%)\n";
  }
  OPENMVG_LOG_INFO << "=================================================================";

  return htmlFileStream.good();
}

} // namespace sfm
} // namespace openMVG
