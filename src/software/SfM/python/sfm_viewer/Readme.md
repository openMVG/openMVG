# ReRun sfm_data.json visualizer sample

## Requirements

`pip install rerun-sdk`

## Demo sample

```
cd src/software/SfM/python/sfm_viewer/sample_data
git clone https://github.com/openMVG/ImageDataset_SceauxCastle.git
cp -r ./ImageDataset_SceauxCastle/images ./

python3 ../rerun_viewer.py --sfm_data ./sfm_data.json --show_keypoints 100_7100.JPG 100_7101.JPG --show_images 100_7100.JPG 100_7106.JPG

```
It will display all the cameras poses, the 3D structures and 2 images, and the keypoints used for the reconstruction on 2 given images.

> **_NOTE:_** Given your needs, you can combine the option `show_keypoints, show_images, show_rays`.

ProTip:
- An openMVG sfm_data.bin can be converted to sfm_data.json by calling ``./openMVG_main_ConvertSfM_DataFormat -i 'sfm_data_in.bin' -o 'sfm_data.json'``

### Display camera to scene rays

`python3 ../rerun_viewer.py --sfm_data ./sfm_data.json --show_rays 100_7103.JPG`
