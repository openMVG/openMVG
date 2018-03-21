**********************************
openMVG_main_SfMInit_ImageListingFromKnownPoses
**********************************

Load datasets with ground truth position and orientation information and pack them into sfm_data.json.
Now support Strcha's Dataset, MiddleBury's Dataset, DTU MVS Dataset, ETH 3D Dataset and Kitti Odometry Dataset.

Note
========================
1. openMVG_main_SfMInit_ImageListingFromKnownPoses only support recitified images nowadays. So when you use the  "DTU MVS Dataset" and "ETH 3D Dataset", you have to pay attention to download the correct images.
2. Please put the ground truth information and images in separate folder so that the  openMVG_main_SfMInit_ImageListingFromKnownPoses can load the correct ground truth file and images.


.. code-block:: c++

  $ openMVG_main_SfMInit_ImageListing -i [] -g [] -t []

Arguments description:

**Required parameters:**

  - **[-i|--imageDirectory]**

    - Path store the images

  - **[-g|--groundTruthDirectory]** 

    - Path store the ground truth(must different with the imageDirectory )

  - **[-t|--groundTruthDataset]**

    - Type of dataset

      -1: Strcha's Dataset
	  -2: MiddleBury's Dataset
      -3: DTU MVS Dataset
	  -4: ETH 3D Dataset
	  -5: Kitti Odometry Dataset

  - **[-o|--outputDirectory]**
    - Path store the sfm_data.json

Demo: Load From Strecha's Dataset
========================
1. Download one of the Strecha's datasets from https://icwww.epfl.ch/~marquez/multiview/denseMVS.html

2. Put the images and ground truth file in separate folder

3. Launch the openMVG_main_SfMInit_ImageListingFromKnownPoses 

.. code-block:: c++

  $ openMVG_main_SfMInit_ImageListingFromKnownPoses -i /home/user/Dataset/Strcha/images -g /home/user/Dataset/Strcha/gt -t 1 -o /home/user/Dataset/Strcha/result
  
4. A sfm_data.json file will be produced that is used by openMVG as a scene description.




