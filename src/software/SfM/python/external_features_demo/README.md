**About**

The Kornia demo script will extract DISK features and match them with LightGlue as well as export features, descriptors, and matches.

**Install dependencies**
```
$ pip install -r requirements.txt
```

**How to use**
```
$ python kornia_demo.py --input <string> --matches <string> --output <string>
```
To see all options:
```
$ python kornia_demo.py -h
```

**Example**

Dataset: https://github.com/openMVG/ImageDataset_SceauxCastle
```
$ openMVG_main_SfMInit_ImageListing -i "...\ImageDataset_SceauxCastle\images" -o "...\ImageDataset_SceauxCastle\images\sfm\matches" -d "...\sensor_width_camera_database.txt"
$ python kornia_demo.py --input "...\ImageDataset_SceauxCastle\images\sfm\matches\sfm_data.json" --matches "...\ImageDataset_SceauxCastle\images\sfm\matches" --output "...\ImageDataset_SceauxCastle\images\sfm\matches\matches.putative.bin"
```
Afterwards, run openMVG_main_GeometricFilter and openMVG_main_SfM as normal.
