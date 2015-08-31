*****************************************
OpenMVS Open Multiple View Stereovision
*****************************************

`OpenMVS <http://cdcseacave.github.io/openMVS/>`_ allows to compute dense points cloud, surface and textured surfaces of OpenMVG scenes.
OpenMVS uses OpenMVG scene thanks to a scene importer.

**Importation**

.. code-block:: c++

  # Import the OpenMVG scene to the OpenMVS data format
  $ InterfaceOpenMVG2 -i PATH/sfm_data.json -o scene.mvs

**Dense point-cloud reconstruction** for obtaining a complete and accurate as possible point-cloud

.. code-block:: c++

  # Compute dense depth map per view and merge the depth map into a consistent point cloud
  $ DensifyPointCloud scene.mvs

**Mesh reconstruction** for estimating a mesh surface that explains the best the input point-cloud

.. code-block:: c++

  # The initial point cloud be:
  # - the calibration one (scene.mvs),
  $ ReconstructMesh scene.mvs
  # - or the dense one (scene_dense.mvs)
  $ ReconstructMesh scene_dense.mvs

**Mesh texturing** for computing a texture to color the mesh

.. code-block:: c++

  # Compute the texture
  $ TextureMesh scene_dense_mesh.mvs
