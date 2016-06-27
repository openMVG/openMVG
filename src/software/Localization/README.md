## Camera Localization

`openMVG_main_cameraLocalization` is used to localize a camera. It can take as 
input a video, a list of images in the form of a simple txt file or a json file. 
It can localize a camera using natural features (SIFT), circular markers (CCTAG) or
both of them (SIFT_CCTAG) [depending on you build options].

### Usage:

```
   -h [ --help ]                      Print this message
  --descriptors arg (=SIFT)          Type of descriptors to use 
                                     {SIFT,CCTAG,SIFT_CCTAG}
  
  --preset arg (=NORMAL)             Preset for the feature extractor when 
                                     localizing a new image 
                                     {LOW,MEDIUM,NORMAL,HIGH,ULTRA}
  --calibration arg                  Calibration file
  --sfmdata arg                      The sfm_data.json kind of file generated 
                                     by OpenMVG.
  --descriptorPath arg               Folder containing the .desc.
  --mediafile arg                    The folder path or the filename for the 
                                     media to track
  --refineIntrinsics                 Enable/Disable camera intrinsics 
                                     refinement for each localized image
  --reprojectionError arg (=4)       Maximum reprojection error (in pixels) 
                                     allowed for resectioning. If set to 0 it 
                                     lets the ACRansac select an optimal value.

  --nbImageMatch arg (=4)            [voctree] Number of images to retrieve in 
                                     database
  --maxResults arg (=10)             [voctree] For algorithm AllResults, it 
                                     stops the image matching when this number 
                                     of matched images is reached. If 0 it is 
                                     ignored.
  --commonviews arg (=3)             [voctree] Number of minimum images in 
                                     which a point must be seen to be used in 
                                     cluster tracking
  --voctree arg                      [voctree] Filename for the vocabulary tree
  --voctreeWeights arg               [voctree] Filename for the vocabulary tree
                                     weights
  --algorithm arg (=AllResults)      [voctree] Algorithm type: FirstBest, 
                                     BestResult, AllResults, Cluster

  --nNearestKeyFrames                [cctag] Number of images to retrieve in the 
                                     database

  --globalBundle                     [bundle adjustment] If --refineIntrinsics 
                                     is not set, this option allows to run a 
                                     final global budndle adjustment to refine 
                                     the scene
  --noDistortion                     [bundle adjustment] It does not take into 
                                     account distortion during the BA, it 
                                     consider the distortion coefficients all 
                                     equal to 0
  --noBArefineIntrinsics             [bundle adjustment] It does not refine 
                                     intrinsics during BA
  --minPointVisibility arg (=0)      [bundle adjustment] Minimum number of 
                                     observation that a point must have in 
                                     order to be considered for bundle 
                                     adjustment
  --visualDebug arg                  If a directory is provided it enables 
                                     visual debug and saves all the debugging 
                                     info in that directory
  --output arg (=trackedcameras.abc) Filename for the SfM_Data export file 
                                     (where camera poses will be stored). 
                                     Default : trackedcameras.abc. It will also
                                     save the localization results (raw data) 
                                     as .json with the same name
```

The `calibration` option expects a plain text file containing the internal camera 
parameters in the following format (as output by `main_cameraCalibration.cpp`):
```
 int #image width
 int #image height
 double #focal
 double #ppx principal point x-coord
 double #ppy principal point y-coord
 double #k0 first distortion parameter
 double #k1
 double #k2
```
