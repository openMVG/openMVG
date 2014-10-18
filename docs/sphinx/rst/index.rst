.. openMVG documentation master file, created by
   sphinx-quickstart on Wed Oct 30 11:05:58 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

======================
openMVG documentation!
======================

.. toctree::
   :maxdepth: 1
   :hidden:
   
   openMVG/openMVG.rst
   openMVG_Samples/openMVG_Samples.rst
   software/softwares.rst
   nonFree/patented.rst
   dependencies/external_libs.rst
   third_party/third_party.rst
   FAQ/FAQ.rst
   bibliography


Introduction
--------------

OpenMVG (Multiple View Geometry) is a library for computer-vision scientists and especially
targeted to the Multiple View Geometry community. It is designed to provide an easy access to the
classical problem solvers in Multiple View Geometry and solve them accurately.


Why another library
--------------

The openMVG credo is **"Keep it simple, keep it maintainable".**
OpenMVG targets readable code that is easy to use and modify by the community.

All the features and modules are unit tested. This test driven development ensures that the code
works as it should and enables more consistent repeatability. Furthermore, it makes it easier for the
user to understand and learn the given features.

openMVG library overview
-------------------------

The openMVG library is cut in various modules:

* **Libraries, core**,

  * comes with unit tests that assert algorithms results and show how use the code.
  
* **Samples**,

  * show how to use the library to build high_level algorithms.
  
* **Binaries**,

  * softwares build to perform toolchain processing

    * features matching in un-ordered photo collection,
    * SfM: Structure from Motion,
    * color harmonization of photo collection.


Acknowledgements
-----------------
openMVG authors would like to thank:

- libmv authors for providing an inspiring base to design the openMVG library.
- Mikros Image and LIGM-Imagine laboratory for support and authorization to make this library as an open-source project.

License
--------------
openMVG library is release under the MPL2 (Mozilla Public License 2.0). It integrates some sub-part
under the MIT (Massachusetts Institute of Technology) and the BSD (Berkeley Software Distribution) license.
Please refer to the license file contained in the source for complete license description.


Dependencies
--------------
OpenMVG come as a standalone distribution, you don't need to install libraries to make it compiles
and run.
On Linux the library will use if available the local png, zlib and jpeg libraries.
   
Indices and tables
----------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

