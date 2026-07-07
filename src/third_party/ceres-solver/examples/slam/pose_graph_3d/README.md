Pose Graph 3D
----------------

The Simultaneous Localization and Mapping (SLAM) problem consists of building a
map of an unknown environment while simultaneously localizing against this
map. The main difficulty of this problem stems from not having any additional
external aiding information such as GPS. SLAM has been considered one of the
fundamental challenges of robotics. A pose graph optimization problem is one
example of a SLAM problem.

The example also illustrates how to use Eigen's geometry module with Ceres'
automatic differentiation functionality. To represent the orientation, we will
use Eigen's quaternion which uses the Hamiltonian convention but has different
element ordering as compared with Ceres's rotation representation. Specifically
they differ by whether the scalar component q_w is first or last; the element
order for Ceres's quaternion is [q_w, q_x, q_y, q_z] where as Eigen's quaternion
is [q_x, q_y, q_z, q_w].

This package defines the necessary Ceres cost functions needed to model the
3-dimensional pose graph optimization problem as well as a binary to build and
solve the problem. The cost functions are shown for instruction purposes and can
be speed up by using analytical derivatives which take longer to implement.


Running
-----------
This package includes an executable `pose_graph_3d` that will read a problem
definition file. This executable can work with any 3D problem definition that
uses the g2o format with quaternions used for the orientation representation. It
would be relatively straightforward to implement a new reader for a different
format such as TORO or others. `pose_graph_3d` will print the Ceres solver full
summary and then output to disk the original and optimized poses
(`poses_original.txt` and `poses_optimized.txt`, respectively) of the robot in
the following format:
```
pose_id x y z q_x q_y q_z q_w
pose_id x y z q_x q_y q_z q_w
pose_id x y z q_x q_y q_z q_w
...
```
where `pose_id` is the corresponding integer ID from the file definition. Note,
the file will be sorted in ascending order for the `pose_id`.

The executable `pose_graph_3d` has one flag `--input` which is the path to the
problem definition. To run the executable,
```
/path/to/bin/pose_graph_3d --input /path/to/dataset/dataset.g2o
```

A script is provided to visualize the resulting output files. There is also an
option to enable equal axes using ```--axes_equal```.
```
/path/to/repo/examples/slam/pose_graph_3d/plot_results.py --optimized_poses ./poses_optimized.txt --initial_poses ./poses_original.txt
```
