**********************************
openMVG_main_SfMInit_ImageListing
**********************************

The first task to process an image dataset in openMVG pipelines consist in creating a **sfm_data.json** file that describe the used image dataset.

This structured file lists for each images an object called a View.
This view store image information and lists:

- the image name,
- the image size,
- the internal camera calibration information (intrinsic parameters) (if any).

Each **View** is associated to a camera **Pose** index and an **Intrinsic** camera parameter group (if any). Group means that camera intrinsic parameters can be shared between some Views (that leads to more stable parameters estimation).

.. code-block:: c++

  $ openMVG_main_SfMInit_ImageListing -i [] -d [] -o []

Arguments description:

**Required parameters:**

  - **[-i|--imageDirectory]**

  - **[-d|--sensorWidthDatabase]** openMVG/src/openMVG/exif/sensor_width_database/sensor_width_camera_database.txt

  - **[-o|--outputDirectory]**

**Optional parameters:**

  - **[-f|--focal]** (value in pixels)

  - **[-k|--intrinsics]** Kmatrix: \"f;0;ppx;0;f;ppy;0;0;1\"

  - **[-c|--camera_model]** Camera model type:

    - 1: Pinhole
    - 2: Pinhole radial 1
    - 3: Pinhole radial 3 (default)
  
  - **[-g|--group_camera_model]**

    - 0-> each view have it's own camera intrinsic parameters
    - 1-> (default) view can share some camera intrinsic parameters

.. code-block:: c++

  // Example
  $ openMVG_main_SfMInit_ImageListing -d /home/user/Dev/openMVG/src/openMVG/exif/sensor_width_database/sensor_width_camera_database.txt -i /home/user/Dataset/ImageDataset_SceauxCastle/images -o /home/user/Dataset/ImageDataset_SceauxCastle/matches
  
  If you have installed OpenMVG on your machine your could do:
  $ openMVG_main_SfMInit_ImageListing -d /usr/share/openMVG/sensor_width_camera_database.txt -i /home/user/Dataset/ImageDataset_SceauxCastle/images -o /home/user/Dataset/ImageDataset_SceauxCastle/matches

It will produce you a sfm_data.json file that is used by openMVG as a scene description.

Once your have computed your dataset description you can compute the image features:

.. toctree::
   :maxdepth: 1

   ./ComputeFeatures.rst

Automatic pixel focal computation from JPEG EXIF metadata
----------------------------------------------------------

Intrinsic image analysis is made from JPEG EXIF metadata and a database of camera sensor width.

**Note:** If you have image with no metadata you can specify the known pixel focal length value directly or let the SfM process find it automatically (at least two images that share common keypoints and with valid intrinsic group must be defined)

The focal in pixel is computed as follow (if the EXIF camera model and maker is found in the provieded sensor width database )

.. math::
      
  \text{focal}_{pix} = \frac{max( w_\text{pix}, h_\text{pix} ) * \text{focal}_\text{mm}} {\text{ccdw}_\text{mm}}

- :math:`\text{focal}_{pix}` the EXIF focal length (pixels),
- :math:`\text{focal}_{mm}` the EXIF focal length (mm),
- :math:`w_\text{pix}, h_\text{pix}` the image of width and height (pixels),
- :math:`\text{ccdw}_\text{mm}` the known sensor width size (mm)

From lists.txt to sfm_data.json
---------------------------------

Old openMVG version (<0.8) use a lists.txt file to describer image parameters.

Example of a lists.txt file where focal is known in advance

.. code-block:: c++

  0000.png;3072;2048;2760;0;1536;0;2760;1024;0;0;1
  0001.png;3072;2048;2760;0;1536;0;2760;1024;0;0;1
  0002.png;3072;2048;2760;0;1536;0;2760;1024;0;0;1
  ...

Example of a lists.txt file where focal is computed from exif data

.. code-block:: c++

  100_7100.JPG;2832;2128;2881.25;EASTMAN KODAK COMPANY;KODAK Z612 ZOOM DIGITAL CAMERA
  100_7101.JPG;2832;2128;2881.25;EASTMAN KODAK COMPANY;KODAK Z612 ZOOM DIGITAL CAMERA
  100_7102.JPG;2832;2128;2881.25;EASTMAN KODAK COMPANY;KODAK Z612 ZOOM DIGITAL CAMERA
  ...

You can convert this file to a valid sfm_data.json file by using the **openMVG_main_ConvertList** application.


