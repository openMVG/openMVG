
*************************************
SfM_data format
*************************************

SfM_Data is a data container. It contains:

- **Views**
  - images
- **Intrinsics**
  - intrinsics camera parameters
- **Poses**
  - extrinsic camera parameters
- **Landmarks**
  - 3D points and their 2D Observations

The View and Intrinsic concept are abstract and so:

- any internal camera model can be used
- any metadata can be stored related to a view.

Dynamic loading of stored object is performed thanks to the cereal serialization library (this library allow polymorphism serialization).

