
********************************************
openMVG_main_ComputeSfM_DataColor
********************************************

Compute the color of the Structure of a sfm_data scene.

Use a very simple approach:

.. code-block:: c++

  a. list the track id with no color
  b. list the most viewed view id
  c. color the track that see the view
  d. go to a. until uncolored track are remaining

Information and usage
========================

The application is designed to run on a sfm_data.json file
The sfm_data file should contains:

- valid view with some defined intrinsics and camera poses,
- (optional existing structure).

  .. code-block:: c++
  
    $ openMVG_main_ComputeSfM_DataColor -i Dataset/out_Reconstruction/sfm_data.json -o Dataset/out_Reconstruction/sfm_data_color.ply

Arguments description:

**Required parameters:**

  - **[-i|--input_file]**

    - a SfM_Data file

  - **[-o|--output_file]**

    - output scene with updated landmarks color

