Export to  OpenMVS Open Multiple View Stereovision
**************************************************

`OpenMVS <http://cdcseacave.github.io/openMVS/>`_ allows to compute dense points cloud, surface and textured surfaces of OpenMVG scenes.

To use OpenMVG scene in OpenMVS you can:
  - use the openMVG_main_openMVG2openMVS exporter interface
  - or use the OpenMVS InterfaceOpenMVG SceneImporter

**Importation/Exportation**

.. code-block:: c++

  # Export the OpenMVG scene to the OpenMVS data format
  $ openMVG_main_openMVG2openMVS -i PATH/sfm_data.bin -d OUTPUT_PATH -o OUTPUT_PATH/Scene

  # Or you can Import the OpenMVG scene to the OpenMVS data format using the OpenMVS binary
  $ InterfaceOpenMVG -i PATH/sfm_data.json -o scene.mvs

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
