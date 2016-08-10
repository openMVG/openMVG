// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm.hpp"
#include "openMVG/exif/exif_IO_EasyExif.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#ifdef HAVE_BOOST
#include <boost/system/error_code.hpp>
#include <boost/filesystem.hpp>
#endif

#include <string>
#include <vector>

using namespace openMVG;
using namespace openMVG::sfm;

/**
 * @brief Update all viewID referenced in the observation of each landmark according 
 * to the provided mapping.
 * 
 * @param[in,out] landmarks The landmarks to update.
 * @param[in] oldIdToNew The mapping between the old ID and the reconmputed UID.
 */
void updateStructureWithNewUID(Landmarks &landmarks, const std::map<std::size_t, std::size_t> &oldIdToNew)
{
  // update the id in the visibility of each 3D point
  for(auto &iter : landmarks)
  {
    Landmark& currentLandmark = iter.second;
    
    // the new observations where to copy the existing ones
    // (needed as the key of the map is the idview)
    Observations newObservations;
    
    for(const auto &iterObs : currentLandmark.obs)
    {
      const auto idview = iterObs.first;
      const Observation &obs = iterObs.second;

      newObservations.emplace(oldIdToNew.at(idview), obs);
    }
    
    assert(currentLandmark.obs.size() == newObservations.size());
    currentLandmark.obs.swap(newObservations);
  }  
}

void regenerateUID(sfm::SfM_Data &sfmdata, std::map<std::size_t, std::size_t> &oldIdToNew, bool sanityCheck = false)
{
  // if the views are empty, nothing to be done. 
  if(sfmdata.GetViews().empty())
    return;
  
  Views newViews;

  for(auto const &iter : sfmdata.views)
  {
    const View& currentView = *iter.second.get();
    const auto &imageName = currentView.s_Img_path;
    
    exif::Exif_IO_EasyExif exifReader(imageName);

    // compute the view UID
    const std::size_t uid = exif::computeUID(exifReader, imageName);

    // update the mapping
    assert(oldIdToNew.count(currentView.id_view) == 0);
    oldIdToNew.emplace(currentView.id_view, uid);
    
    // add the view to the new map using the uid as key and change the id
    assert(newViews.count(uid) == 0);
    newViews.emplace(uid, iter.second);
    newViews[uid]->id_view = uid;
  }
  
  assert(newViews.size() == sfmdata.GetViews().size());
  sfmdata.views.swap(newViews);
  
  
  // update the id in the visibility of each 3D point
  updateStructureWithNewUID(sfmdata.structure, oldIdToNew);
  
   // update the id in the visibility of each 3D point
  updateStructureWithNewUID(sfmdata.control_points, oldIdToNew);
  
  if(!sanityCheck)
    return;
  
  // sanity check
  for(auto &iter : sfmdata.structure)
  {
    Landmark& currentLandmark = iter.second;
    for(const auto &iterObs : currentLandmark.obs)
    {
      const auto idview = iterObs.first;
      const Observation &obs = iterObs.second;

      // there must be a view with that id (in the map) and the view must have 
      // the same id (the member)
      assert(sfmdata.views.count(idview) == 1);
      assert(sfmdata.views.at(idview)->id_view == idview);
    }
  }  
  
  // sanity check
  for(auto &iter : sfmdata.control_points)
  {
    Landmark& currentLandmark = iter.second;
    for(const auto &iterObs : currentLandmark.obs)
    {
      const auto idview = iterObs.first;
      const Observation &obs = iterObs.second;

      // there must be a view with that id (in the map) and the view must have 
      // the same id (the member)
      assert(sfmdata.views.count(idview) == 1);
      assert(sfmdata.views.at(idview)->id_view == idview);
    }
  }
  
}

// Convert from a SfM_Data format to another
int main(int argc, char **argv)
{
  CmdLine cmd;

  std::string sSfM_Data_Filename_In;
  std::string sSfM_Data_Filename_Out;
#ifdef HAVE_BOOST
  std::string matchDir;
#endif

  cmd.add(make_option('i', sSfM_Data_Filename_In, "input_file"));
  cmd.add(make_switch('V', "VIEWS"));
  cmd.add(make_switch('I', "INTRINSICS"));
  cmd.add(make_switch('E', "EXTRINSICS"));
  cmd.add(make_switch('S', "STRUCTURE"));
  cmd.add(make_switch('O', "OBSERVATIONS"));
  cmd.add(make_switch('C', "CONTROL_POINTS"));
  cmd.add(make_switch('u', "uid"));
  cmd.add(make_option('o', sSfM_Data_Filename_Out, "output_file"));
#ifdef HAVE_BOOST
  cmd.add(make_option('m', matchDir, "matchDirectory"));
#endif

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch(const std::string& s) {
      std::cerr << "Usage: " << argv[0] << '\n'
        << "[-i|--input_file] path to the input SfM_Data scene\n"
        << "[-o|--output_file] path to the output SfM_Data scene\n"
        << "\t .json, .bin, .xml, .ply, .baf"
#if HAVE_ALEMBIC
           ", .abc"
#endif
           "\n"
        << "\n[Options to export partial data (by default all data are exported)]\n"
        << "\nUsable for json/bin/xml format\n"
        << "[-V|--VIEWS] export views\n"
        << "[-I|--INTRINSICS] export intrinsics\n"
        << "[-E|--EXTRINSICS] export extrinsics (view poses)\n"
        << "[-S|--STRUCTURE] export structure\n"
        << "[-O|--OBSERVATIONS] export 2D observations associated with 3D structure\n"
        << "[-C|--CONTROL_POINTS] export control points\n"
        << "[-u|--uid] (re-)compute the unique ID (UID) for the views\n"
#ifdef HAVE_BOOST              
        << "[-m|--matchDirectory] the directory containing the features used for the\n"
           "    reconstruction. If provided along the -u option, it creates symbolic\n"
           "    links to the .desc and .feat with the new UID as name. This can be\n"
           "    for legacy reconstructions that were not made using UID"
#endif
        << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  if (sSfM_Data_Filename_In.empty() || sSfM_Data_Filename_Out.empty())
  {
    std::cerr << "Invalid input or output filename." << std::endl;
    return EXIT_FAILURE;
  }

  // OptionSwitch is cloned in cmd.add(),
  // so we must use cmd.used() instead of testing OptionSwitch.used
  int flags =
    (cmd.used('V') ? VIEWS      : 0)
  | (cmd.used('I') ? INTRINSICS : 0)
  | (cmd.used('E') ? EXTRINSICS : 0)
  | (cmd.used('O') ? OBSERVATIONS : 0)
  | (cmd.used('S') ? STRUCTURE  : 0);

  flags = (flags) ? flags : ALL;
  
  const bool recomputeUID = cmd.used('u');

  // Load input SfM_Data scene
  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename_In, ESfM_Data(ALL)))
  {
    std::cerr << std::endl
      << "The input SfM_Data file \"" << sSfM_Data_Filename_In << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }
  
  if(recomputeUID)
  {
    std::cout << "Recomputing the UID of the views..." << std::endl;
    std::map<std::size_t, std::size_t> oldIdToNew;
    regenerateUID(sfm_data, oldIdToNew);
    
#ifdef HAVE_BOOST
    if(!matchDir.empty())
    {
      std::cout << "Generating alias for .feat and .desc with the UIDs" << std::endl;
      for(const auto& iter : oldIdToNew)
      {
        const auto oldID = iter.first;
        const auto newID = iter.second;
        const auto oldFeatfilename = stlplus::create_filespec(matchDir, std::to_string(oldID), ".feat");
        const auto newFeatfilename = stlplus::create_filespec(matchDir, std::to_string(newID), ".feat");
        const auto oldDescfilename = stlplus::create_filespec(matchDir, std::to_string(oldID), ".desc");
        const auto newDescfilename = stlplus::create_filespec(matchDir, std::to_string(newID), ".desc");

        if(!(stlplus::is_file(oldFeatfilename) && stlplus::is_file(oldDescfilename)))
        {
          std::cerr << "Cannot find the features file for view ID " << oldID
                      << std::endl;
          return EXIT_FAILURE;
        }
        boost::system::error_code ec;
        boost::filesystem::create_symlink(oldFeatfilename, newFeatfilename, ec);
        if(ec)
        {
          std::cerr << "Error while creating " << newFeatfilename << ": " << ec.message() << std::endl;
          return EXIT_FAILURE;
        }
        boost::filesystem::create_symlink(oldDescfilename, newDescfilename, ec);
        if(ec)
        {
          std::cerr << "Error while creating " << newDescfilename << ": " << ec.message() << std::endl;
          return EXIT_FAILURE;
        }
      }
    }
#endif
    
  }

  // Export the SfM_Data scene in the expected format
  if (!Save(sfm_data, sSfM_Data_Filename_Out, ESfM_Data(flags)))
  {
    std::cerr << std::endl
      << "An error occured while trying to save \"" << sSfM_Data_Filename_Out << "\"." << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
