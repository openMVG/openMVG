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

As OpenMVG uses some C++11 features, you must have a C++11 ready compiler:

- Visual Studio >= 2015 (recommended)
- GCC >= 4.8.1
- Clang >= 3.3

General information for OpenMVG CMake options:

- OpenMVG_BUILD_TESTS (ON/OFF(default))
    - Build OpenMVG unit tests
- OpenMVG_BUILD_EXAMPLES (ON/OFF(default))
    - Build OpenMVG example applications.

Note: options does not affect binaries under 'software'


General information for OpenMVG SfM pipelines:

- OpenMVG can export graphs as Graphviz .dot files and render them as SVG files. If you want this graph visualization feature, please install Graphviz.

Compilation
----------------

- [Linux](#linux)
- [Windows](#windows)
- [MacOS](#macos)


Linux compilation
-----------------
<a name="linux"></a>

Install the required external libraries.
```shell
$ sudo apt-get install libpng-dev libjpeg-dev libtiff-dev libxxf86vm1 libxxf86vm-dev libxi-dev libxrandr-dev
```
If you want see the view graph svg logs, install Graphviz.
```shell
$ sudo apt-get install graphviz
```

Build OpenMVG.
```shell
$ git clone --recursive https://github.com/openMVG/openMVG.git
$ cd openMVG
$ ls
 AUTHORS BUILD  docs  logo  README  src  ...
$ cd ..
$ mkdir openMVG_Build
$ cd openMVG_Build
```
If you want to add unit tests and examples to the build, run
```shell
$ cmake -DCMAKE_BUILD_TYPE=RELEASE -DOpenMVG_BUILD_TESTS=ON -DOpenMVG_BUILD_EXAMPLES=ON . ../openMVG/src/
```
otherwise
```shell
$ cmake -DCMAKE_BUILD_TYPE=RELEASE . ../openMVG/src/
```
If you want to have an IDE openable project with Code::Blocks:
```shell
$ cmake -G "CodeBlocks - Unix Makefiles" -DCMAKE_BUILD_TYPE=RELEASE . ../openMVG/src/
```

Compile the project
```shell
$ make
```
or for a multi-core compilation. (Replace NBcore with the number of threads)
```shell
$ make -j NBcore
```

Run tests (if requested at CMake step)
```shell
$ make test
```

Have fun with the samples
```shell
$ cd openMVG_Samples
```

Windows compilation
-------------------
<a name="windows"></a>

Checkout the project
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

Mac compilation
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
If you want to add unit tests and examples to the build, run
```shell
$ cmake -DCMAKE_BUILD_TYPE=RELEASE -DOpenMVG_BUILD_TESTS=ON -DOpenMVG_BUILD_EXAMPLES=ON -G "Xcode" . ../openMVG/src/
```

otherwise you can use standard makefiles
```shell
$ cmake -DCMAKE_BUILD_TYPE=RELEASE . ../openMVG/src/
```

Using OpenCV sample
--------------------

Add `-DOpenMVG_USE_OPENCV=ON` to your cmake command and set the OpenCV_DIR variable to your OpenCV build directory
e.g. `-DOpenCV_DIR="/home/user/Dev/github/itseez/opencv_Build" -DOpenMVG_USE_OPENCV=ON`

Using OpenMVG as a third party library dependency in CMake
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

Specify in CMake where OpenMVG has been installed by using the CMake OpenMVG_DIR variable
e.g. `-DOpenMVG_DIR:STRING="YourInstallPath"/share/openMVG/cmake`

A message will be displayed if OpenMVG is found or not at the CMake configure step.

