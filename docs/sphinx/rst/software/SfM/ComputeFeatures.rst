
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
      - ULTRA: !!Can be time consuming!!


**Use mask to filter keypoints/regions**

  Sometime you may want to compute features/regions only on some parts of your images. It could include the following cases: 

  - You know that in your acquisition some areas may disturb the SfM pipeline (even if OpenMVG is known to be robust to a fair amount of outliers) 
    and lower the quality of the subsequent reconstruction. For example, in some close range configurations, you may prefer to move the object itself 
    instead of moving the camera. Masks can also help you to deal with hemispherical fish-eyes, by masking useless zone of the sensor.
  - You want to speed up the computation by reducing the number features/regions and thus the number of tie points.

  For this kind of needs you can use a mask. A mask is simply a binary image having the same size (width and height) than the target image.
  The black areas on a mask denote the "bad parts", *i.e.* the areas to be masked and for which descriptors are not computed. A point is kept if the mask value at the point position is different than 0.
  In openMVG_main_ComputeFeatures, the association of a mask and an image is implicit. It uses the following conventions:

  - It tries to load a global mask.png file from directory where the the SfM container file (sfm_data.*) is stored.
  - It tries to load an individual <image_name>_mask.png from directory where the current image is stored 

  The individual mask **always** takes precedence over the global one 

Once openMVG_main_ComputeFeatures is done you can compute the Matches between the computed description.

.. toctree::
   :maxdepth: 1

   ./ComputeMatches.rst

Export detected regions as SVG files:
----------------------------------------

* **Detected keypoints**: openMVG_main_exportKeypoints


