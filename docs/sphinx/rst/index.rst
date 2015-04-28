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


OpenMVG (Multiple View Geometry) is a library for computer-vision scientists and targeted for the Multiple View Geometry community.

It is designed to provide an easy access to:

* Accurate Multiple View Geometry problem solvers,
* Tiny libraries to perform tasks from feature detection/matching to Structure from Motion,
* Complete Structure from Motion pipelines.


Why another library
-------------------

The openMVG credo is **"Keep it simple, keep it maintainable".**.
It means provide a readable code that is easy to use and modify by the community.

All the **features and modules are unit-tested**.
This test driven development ensures:

* more consistent repeatability (assert code works as it should in time),
* that the code is used in a real context to show how it should be used.

openMVG overview
-------------------------

OpenMVG is cut in various modules:

* **Libraries**, core modules,

  * comes with unit tests that assert algorithms results and show how use the code.
  
* **Samples**,

  * show how to use the library to build high_level algorithms.
  
* **Softwares**,

  *  ready to use tools to perform toolchain processing:

    * features matching in un-ordered photo collection,
    * SfM: tools and Structure from Motion pipelines,
    * color harmonization of photo collection.

Cite Us
=======

If you use openMVG for a publication, please cite it as::

    @misc{openMVG,
      author = "Pierre Moulon and Pascal Monasse and Renaud Marlet and Others",
      title = "OpenMVG",
      howpublished = "\url{https://github.com/openMVG/openMVG}",
    }


Acknowledgements
-----------------
openMVG authors would like to thank:

- libmv authors for providing an inspiring base to design the openMVG library,
- Mikros Image, LIGM-Imagine laboratory and Foxel SA for support,
- Mikros Image, LIGM-Imagine laboratory authorization to make this library an open-source project.

License
--------------
openMVG library is release under the MPL2 (Mozilla Public License 2.0). It integrates some sub-part
under the MIT (Massachusetts Institute of Technology) and the BSD (Berkeley Software Distribution) license.
Please refer to the copyright.md and license files contained in the source for complete license description.


Dependencies
--------------
OpenMVG comes as a standalone distribution, you don't need to install libraries to make it compiles
and run.
On Linux openMVG will use the local png, zlib and jpeg libraries if they are availables.
   
Indices and tables
----------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

