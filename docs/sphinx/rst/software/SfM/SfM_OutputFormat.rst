
*************************************
SfM_Output openMVG format
*************************************

OpenMVG export the SfM data in the following dicrectory tree structure:

.. code-block:: c++

  views.txt
  images [Directory]
  cameras [Directory]
  cameras_disto [Directory]
  clouds [Directory]

-----------------
view.txt
-----------------

A txt file that contains information about the sucessfully oriented cameras:

.. code-block:: c++

  images
  cameras
  3
  100_7100.JPG 2832 2128 100_7100.bin 0.453315 2.82001
  100_7101.JPG 2832 2128 100_7101.bin 0.46713 2.995
  100_7102.JPG 2832 2128 100_7102.bin 0.517922 3.12022

It provides for each camera:

ImageName, Width, Height, cameraName, NearPlane, FarPlane

------------------
images directory
------------------

Contains successfully oriented images as undistorted images.

------------------
cameras directory
------------------

Each P[4x3] camera matrices saved as file (12 double in binary format). One file per camera.
See pmoulon/openMVG/src/openMVG/cameras/Camera_IO.hpp to read it.

.. code-block:: c++

  PinholeCamera cam;
  bool bOk = load("filename.bin", cam);

----------------------------
cameras_disto directory
----------------------------

Each camera saved as BrownPinholeCamera. It contains distortion and camera orientation, position and intrinsics.
Each camera saved in ascii format, one file per camera:

.. code-block:: c++

  focal, ppx, ppy, k1, k2, k3;
  R[3x3]
  t[3x1]
  
  // To read it we suggest use
  BrownPinholeCamera cam;
  bool bOk = load("filename.txt", cam);

------------------
clouds directory
------------------

Contain data relative to the Structure and visibility of the scene in two different file.

**calib.ply** file that can be opened in Meshlab or CloudCompare.
Color field have been updated only on request, else it will be white.

.. code-block:: c++

   X Y Z red green blue confidence #Visibility [ImageIndex]
   float float float uchar uchar uchar float int intList

So for each 3D point of the structure the color and camera visibility is depicted.
The confidence field is not used by now.
Visibility indexes are related to the Inth camera listed in the views.txt file.

**visibility.txt** file contains more detailed information

.. code-block:: c++

  X Y Z #Visibility [ImageIndex, FeatureIndex]
  double double double int intList

Visibility indexes are related to the Inth camera listed in the views.txt file.
Feature indexes are related to the .feat files exported at matching time.

------------------------------------------------------
How to load all the data as a ready to use project
------------------------------------------------------

See openMVG/software/SfMViewer/document.h

.. code-block:: c++

    std::string openMVG_OutputPath = "Foo/SfM_Output";
    Document m_doc;    
    bool bOk = m_doc.load(openMVG_OutputPath);
    

