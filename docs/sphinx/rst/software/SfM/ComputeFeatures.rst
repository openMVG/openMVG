
***************************************************************
openMVG_main_ComputeFeatures
***************************************************************

Compute image description for a given sfm_data.json file.
For each view it compute the image description (local regions) and store them on disk

.. code-block:: c++

  $ openMVG_main_ComputeFeatures -i [..\matches\sfm_data.json] -o [...\matches]

Arguments description:

**Required parameters:**

  - **[-i|--input_file]**

    - a SfM_Data file

  - **[-o|--outdir path]**
    
    - path were image description will be stored

**Optional parameters:**
     
  - **[-f|--force: Force to recompute data]**

    - 0: (default) reload previously computed data (useful when you have kill the process and want to continue to compute)
    - 1: useful when you change have changed a command line parameter, force recomputing and re-saving.

  - **[-m|--describerMethod]**

    - Used method to describe an image:

      - SIFT: (default),
      - AKAZE_FLOAT: AKAZE with floating point descriptors,
      - AKAZE_MLDB:  AKAZE with binary descriptors.

  - **[-u|--upright]**

    - Use Upright feature or not

      - 0: (default, rotation invariance) 
      - 1: extract upright feature (orientation angle = 0°)

  - **[-p|--describerPreset]**

    - Used to control the Image_describer configuration:

      - NORMAL,
      - HIGH,
      - ULTRA: !!Can be time consumming!!

Once openMVG_main_ComputeFeatures is done you can compute the Matches between the computed description.

.. toctree::
   :maxdepth: 1

   ./ComputeMatches.rst

Export detected regions as SVG files:
----------------------------------------

* **Detected keypoints**: openMVG_main_exportKeypoints


