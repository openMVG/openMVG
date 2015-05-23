
********************************************
openMVG_main_ExportUndistortedImages
********************************************

This application export undistorted images from known camera parameter intrinsic.

Algorithm of the application

.. code-block:: c++

  Require: internal + camera calibration
  Require: images
  Ensure: undistorted images
  for each view
      if view has valid intrinsic
        undistort and save the undistorted view

Information and usage
========================

The chain is designed to run on a sfm_data.json file.
The sfm_data file should contains:
- valid view with some defined intrinsics,

  .. code-block:: c++
  
    $ openMVG_main_ExportUndistortedImages -i Dataset/out_Reconstruction/sfm_data.json -o Dataset/out_Reconstruction/undistortedImages

Arguments description:

**Required parameters:**

  - **[-i|--input_file]**

    - a SfM_Data file with valid intrinsics and poses and optional structure

  - **[-o|--outdir]**

    - path where the undistorted images will be stored


