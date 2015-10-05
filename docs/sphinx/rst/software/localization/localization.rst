***************************
Localization
***************************

This module provide tools to localize an image into an existing Structure from Motion scene.

openMVG_main_SfM_Localization
==============================

.. code-block:: c++

  $ openMVG_main_SfM_Localization -i [] -m [] -o [] -q []

Arguments description:

**Required parameters:**

  - **[-i|--input_file]** The input SfM_Data scene (must contains a structure and camera poses)

  - **[-m|--match_dir]** path to the matches that corresponds to the provided SfM_Data scene

  - **[-o|--out_dir]** path to save the found camera position

**Optional parameters:**

  - **[-q|--query_image]** The query image to locate

.. code-block:: c++

  // Example
  $ openMVG_main_SfM_Localization -i /home/user/Dataset/ImageDataset_SceauxCastle/reconstruction/sfm_data.json -o /home/user/Dataset/ImageDataset_SceauxCastle/matches -o ./ -q /home/user/Dataset/ImageDataset_SceauxCastle/images/100_7100.JPG
