Pose Graph 2D
----------------

The Simultaneous Localization and Mapping (SLAM) problem consists of building a
map of an unknown environment while simultaneously localizing against this
map. The main difficulty of this problem stems from not having any additional
external aiding information such as GPS. SLAM has been considered one of the
fundamental challenges of robotics. A pose graph optimization problem is one
example of a SLAM problem.

This package defines the necessary Ceres cost functions needed to model the
2-dimensional pose graph optimization problem as well as a binary to build and
solve the problem. The cost functions are shown for instruction purposes and can
be speed up by using analytical derivatives which take longer to implement.

Running
-----------
This package includes an executable `pose_graph_2d` that will read a problem
definition file. This executable can work with any 2D problem definition that
uses the g2o format. It would be relatively straightforward to implement a new
reader for a different format such as TORO or others. `pose_graph_2d` will print
the Ceres solver full summary and then output to disk the original and optimized
poses (`poses_original.txt` and `poses_optimized.txt`, respectively) of the
robot in the following format:

```
pose_id x y yaw_radians
pose_id x y yaw_radians
pose_id x y yaw_radians
...
```

where `pose_id` is the corresponding integer ID from the file definition. Note,
the file will be sorted in ascending order for the `pose_id`.

The executable `pose_graph_2d` has one flag `--input` which is the path to the
problem definition. To run the executable,

```
/path/to/bin/pose_graph_2d --input /path/to/dataset/dataset.g2o
```

A python script is provided to visualize the resulting output files.
```
/path/to/repo/examples/slam/pose_graph_2d/plot_results.py --optimized_poses ./poses_optimized.txt --initial_poses ./poses_original.txt
```
