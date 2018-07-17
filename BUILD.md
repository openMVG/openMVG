OpenMVG (Open Multiple View Geometry)
=====================================

Build instructions
------------------

Required tools:

- CMake
- Git
- C/C++ compiler (GCC, Visual Studio or Clang)

Optional tools:

- QT version >= v5.4

As OpenMVG uses some C++11 features, you must have a C++11 ready compiler:

- Visual Studio >= 2015 (recommended)
- GCC >= 4.8.1
- Clang >= 3.3

Note:

- CMAKE OpenMVG variables you can configure:

  - OpenMVG_BUILD_TESTS (ON/OFF(default))
      - Build OpenMVG unit tests
  - OpenMVG_BUILD_EXAMPLES (ON/OFF(default))
      - Build OpenMVG example applications.
  - OpenMVG_BUILD_SOFTWARES (ON(default)/OFF)
      - Build OpenMVG applications.

- General information for OpenMVG SfM pipelines:

  - OpenMVG can export graphs as Graphviz .dot files and render them as SVG files. If you want this graph visualization feature, please install Graphviz.

Checking out the project and build it
--------------------------------------

- [Getting the project](#checkout)
- [Compiling on Linux](#linux)
- [Compiling on Windows](#windows)
- [Compiling on MacOS](#macos)
- [Compiling using VCPKG](#vcpkg)

Getting the project
--------------------
<a name="checkout"></a>

Getting the sources (and the submodules):
```shell
$ git clone --recursive https://github.com/openMVG/openMVG.git
```
or
```shell
$ git clone https://github.com/openMVG/openMVG.git
$ cd openMVG
$ git submodule init
$ git submodule update
```

Compiling on Linux
-------------------
<a name="linux"></a>

1. Install the required external libraries.
```shell
$ sudo apt-get install libpng-dev libjpeg-dev libtiff-dev libxxf86vm1 libxxf86vm-dev libxi-dev libxrandr-dev
```
If you want see the view graph svg logs, install Graphviz.
```shell
$ sudo apt-get install graphviz
```

2. Checkout OpenMVG.
```shell
$ git clone --recursive https://github.com/openMVG/openMVG.git
$ mkdir openMVG_Build && cd openMVG_Build
```

3. Configure and build
```shell
$ cmake -DCMAKE_BUILD_TYPE=RELEASE ../openMVG/src/
$ cmake --build . --target install
```

Run tests using make or ctest (if requested in the CMake command line with `-DOpenMVG_BUILD_TESTS=ON`)
```shell
$ make test
$ ctest --output-on-failure -j
```

Compiling on Windows
---------------------
<a name="windows"></a>

1. Checkout the project
```shell
$ git clone --recursive https://github.com/openMVG/openMVG.git
```

1. Open cmake-gui.
2. Fill the source path with the src OpenMVG path.
3. Fill the build path with a new directory.
4. Select your Visual Studio IDE and click configure and then generate.
5. Open the .sln solution created in your build directory.
6. Change the target to Release.
7. Compile the libraries and binaries samples.

Another options is to build with cmake in console. I recommend launching `VS2015/VS2017 x64 Native Tools Command Prompt` which has build environment already set up or to use the [VCPKG](#vcpkg) alternative.

Compiling on Mac
-------------------
<a name="macos"></a>

```shell
$ git clone --recursive https://github.com/openMVG/openMVG.git
$ cd openMVG
$ ls
 AUTHORS BUILD  docs  logo  README  src  ...
$ cd ..
$ mkdir openMVG_Build
$ cd openMVG_Build
```

If you want to use Xcode and compile by using the command line, run
```
$ cmake -DCMAKE_BUILD_TYPE=RELEASE -G "Xcode" . ../openMVG/src/
$ xcodebuild -configuration Release
```
otherwise you can use standard makefiles
```shell
$ cmake -DCMAKE_BUILD_TYPE=RELEASE . ../openMVG/src/
```

Compiling using VCPKG
-----------------------
<a name="vcpkg"></a>

Checking and build VCPKG
```shell
$ git clone https://github.com/Microsoft/vcpkg
$ cd vcpkg
$ ./bootstrap-vcpkg.sh/bat
```

Checking OpenMVG 3rd party dependencies by using VCPKG (configure your build triplet if needed see `--triplet`)
```shell
$ ./vcpkg install cereal ceres eigen3 libjpeg-turbo libpng tiff 
```
Note: If you want to use ceres with a sparse back end, please use one of this choice `ceres[cxsparse]` or `ceres[suitesparse]` or `ceres[eigensparse]`.

Checking out OpenMVG and build it
```shell
$ git clone --recursive https://github.com/openMVG/openMVG.git
$ mkdir openMVG_Build
$ cd openMVG_Build
$ cmake ../openMVG/src/ -DCMAKE_TOOLCHAIN_FILE=<VCPK_ROOT>/scripts/buildsystems/vcpkg.cmake
$ cmake --build .
```

Using OpenCV sample
--------------------

Add `-DOpenMVG_USE_OPENCV=ON` to your cmake command and set the OpenCV_DIR variable to your OpenCV build directory
e.g. `-DOpenCV_DIR="/home/user/Dev/github/itseez/opencv_Build" -DOpenMVG_USE_OPENCV=ON`

Using OpenMVG as a third party library dependency with CMake
-------------------------------------------------------------

OpenMVG can be used as a third party library once it has been installed.
Because it can use its own Ceres version, it's better to install it locally and not in system files.
So please consider using the CMAKE_INSTALL_PREFIX CMake variable to specify a local installation directory.

Here is the syntax to add the variable to the cmake command (use absolute path):
`-DCMAKE_INSTALL_PREFIX:STRING="YourInstallPath"`
e.g `-DCMAKE_INSTALL_PREFIX:STRING="/home/user/Dev/github/openMVG_Build/openMVG_install"`

Perform `make` and `make install`.

Once the library has been installed, go to your project that wants to use OpenMVG as an external library and add

```
find_package(OpenMVG REQUIRED)
include_directories(${OPENMVG_INCLUDE_DIRS})
add_executable(main main.cpp)
target_link_libraries(main ${OPENMVG_LIBRARIES})
```

or with modern target-based approach (CMake 3.0+, Includes directories will be added by CMake target transitivity)

```
find_package(OpenMVG REQUIRED)
add_executable(main main.cpp)
target_link_libraries(main
  PRIVATE
    OpenMVG::openMVG_sfm
    OpenMVG::openMVG_matching
    ...
)
```

The list of available openMVG libraries is (use only the ones you need):
```
OpenMVG::openMVG_camera
OpenMVG::openMVG_exif
OpenMVG::openMVG_features
OpenMVG::openMVG_geodesy
OpenMVG::openMVG_geometry
OpenMVG::openMVG_graph
OpenMVG::openMVG_image
OpenMVG::openMVG_linearProgramming
OpenMVG::openMVG_matching
OpenMVG::openMVG_matching_image_collection
OpenMVG::openMVG_multiview
OpenMVG::openMVG_numeric
OpenMVG::openMVG_robust_estimation
OpenMVG::openMVG_sfm
OpenMVG::openMVG_system
```

If OpenMVG has been installed by using the CMake OpenMVG_DIR variable you can specify where the install have been done manually by using:

  `-DOpenMVG_DIR:STRING="YourInstallPath"/share/openMVG/cmake`

A message will be displayed if OpenMVG is found or not at the CMake configure step.
