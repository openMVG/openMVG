Export to MVE (Multi-View Environment)
***************************************

`MVE <https://github.com/simonfuhrmann/mve>`_ can import a converted openMVG SfM scene and use it to create dense depth map and complete dense 3D models.

.. code-block:: c++

  # Convert the openMVG SfM scene to the MVE format
  $ openMVG_main_openMVG2MVE -i Dataset/outReconstruction/sfm_data.bin -o Dataset/outReconstruction

  #--
  # shell script example
  #--

  directory=Dataset/outReconstruction/MVE
  resolution=2

  # MVE
  dmrecon -s$resolution $directory
  scene2pset -ddepth-L$resolution -iundist-L$resolution -n -s -c $directory $directory/OUTPUT.ply
  fssrecon $directory/OUTPUT.ply $directory/OUTPUT_MESH.ply
  meshclean $directory/OUTPUT_MESH.ply $directory/OUTPUT_MESH_CLEAN.ply

  Call any tool without arguments to see the help.

Export to MVS Texturing
=======================

If you don't want to use the full MVE pipeline but only `MVS Texturing <https://github.com/nmoehrle/mvs-texturing>`_ [Waechter2014]_ to project a set of oriented images on a mesh, one solution is to use the openMVG_main_openMVG2MVSTEXTURING binary. This binary converts your SfM_Data file into one format used by MVS Texturing. In addition, you may need to undistort your images with openMVG_main_ExportUndistortedImages as it's not handled by the openMVG_main_openMVG2MVSTEXTURING tool.
