*****************************
Intrinsic groups
*****************************

======================================
openMVG intrinsic group explaination
======================================

openMVG groups pictures that share common intrinsic parameter to make parameters estimation more stable.
So each camera will have it's own extrinsic parameter group but can share intrinsic parameters with an image collection.

Intrinsic image analysis is made from JPEG EXIF metadata and a database of camera sensor width.
If you have image with no metadata you can specify the known pixel focal length value directly.

Intrinsic analysis and export:
-----------------------------------

The process exports in outputDirectory/**lists.txt** file the extracted camera intrinsics, one line per image.

  .. code-block:: c++

    $ cd ./software/SfM/Release/
    $ openMVG_main_CreateList [-i|--imageDirectory] [-d|--sensorWidthDatabase] [-o|--outputDirectory] [-f|--focal]

  - Usage of the automatic chain (with JPEG images)
  
  .. code-block:: c++
  
    $ openMVG_main_CreateList -i Dataset/images/ -o Dataset/matches/ -d ./software/SfM/cameraSensorWidth/cameraGenerated.txt


  The focal in pixel is computed according the following formula:

    .. math::
      
      \text{focal}_{pix} = \frac{max( w_\text{pix}, h_\text{pix} ) * \text{focal}_\text{mm}} {\text{ccdw}_\text{mm}}

    - :math:`\text{focal}_{pix}` the EXIF focal length (pixels),
    - :math:`\text{focal}_{mm}` the EXIF focal length (mm),
    - :math:`w_\text{pix}, h_\text{pix}` the image of width and height (pixels),
    - :math:`\text{ccdw}_\text{mm}` the known sensor width size (mm)

      - used if the EXIF camera model and maker is found in the sensorWidthDatabase database openMVG/src/software/SfM/cameraSensorWidth/cameraGenerated.txt

  - If all the camera have the same focal length or no EXIF info, you can set the pixel focal length directly
  
  .. code-block:: c++
  
    $ openMVG_main_CreateList -f 2750(i.e in pixels) -i Dataset/images/ -o Dataset/matches/

Depending if EXIF data and if the camera model and make is registered in the camera databset, per lines the following information are displayed:

  - no information can be extracted (only the principal point position is exported)
  - the image contain EXIF Jpeg approximated focal length
    
    - if the camera sensor is saved in the openMVG database the focal length and principal point is exported
    - else no focal camera can be computed, only the image size is exported.


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
