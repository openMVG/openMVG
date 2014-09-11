=====================================
OpenMVG (open Multiple View Geometry)
=====================================

-----------------
Build instruction
-----------------

Required tools:
* Cmake 
* Git
* c/c++ compiler (gcc or visual studio or clang)

Getting the sources:
$ git clone --recursive https://github.com/openMVG/openMVG.git
or
$ git clone https://github.com/openMVG/openMVG.git
$ cd openMVG
$ git submodule init
$ git submodule update

------------------------------------
Using as library dependency in cmake
------------------------------------
Adding following lines to your CMakeLists.txt should provide OpenMVG usable as
static library:

 add_subdirectory(openMVG/src)
 include_directories(${OpenMVG_INCLUDES})
 target_link_libraries(target ${OpenMVG_LIBS})

Information about required dependencies, standalone build and platform 
specificity can be found below.

--------------------------
General informations
for openMVG cmake options
--------------------------
OpenMVG_BUILD_TESTS (ON/OFF(default))=> Build openMVG unit tests
OpenMVG_BUILD_EXAMPLES (ON/OFF(default))=> Build OpenMVG example applications.
  Does not affect binaries under 'software'

-----------------
Linux compilation
-----------------

Setup the required external library.
* sudo apt-get install libpng-dev libjpeg-dev libtiff-dev libxxf86vm1 libxxf86vm-dev libxi-dev libxrandr-dev

 $ git clone --recursive https://github.com/openMVG/openMVG.git
 $ cd openMVG
 $ ls
  AUTHORS BUILD  docs  logo  README  src  ...
 $ cd ..
 $ mkdir openMVG_Build
 $ cd openMVG_Build
 $ cmake -DCMAKE_BUILD_TYPE=RELEASE . ../openMVG/src/
If you want enable unit tests and examples to the build:
 $ cmake -DCMAKE_BUILD_TYPE=RELEASE -DOpenMVG_BUILD_TESTS=ON -DOpenMVG_BUILD_EXAMPLES=ON . ../openMVG/src/

=> In order to use the MOSEK 6 back-end for the linear programming openMVG module
  - Check that you have an up-to-date MOSEK licence, else openMVG MOSEK unit test will fail.

 $ cmake -DCMAKE_BUILD_TYPE=RELEASE
    -DMOSEK_SEARCH_HEADER="~/Documents/Lib/mosek/6/tools/platform/linux64x86/h"
    -DMOSEK_SEARCH_LIB="~/Documents/Lib/mosek/6/tools/platform/linux64x86/bin"
    . ../openMVG/src/


If you want have an IDE openable project with codeblocks:
 $ cmake -G "CodeBlocks - Unix Makefiles" -DCMAKE_BUILD_TYPE=RELEASE . ../openMVG/src/

Compile the project
 $ make

For a multi-core compilation (Replace NBcore with the number of threads)
 $ make -j NBcore

Launch test (if asked at cmake step)
 $ make test

Have fun with the samples
 $ cd openMVG_Samples

-------------------
Windows compilation
-------------------

Checkout the project
 $ git clone --recursive https://github.com/openMVG/openMVG.git

Open cmake-gui
Fill the source path with the src openMVG path.
Fill the build path with a new directory
Select your Visual Studio IDE and click configure and then generate

Open the .sln solution created in your build directory.
Change the target to Release.
Compile the libraries and binaries samples.

-------------------
Mac compilation
-------------------
 $ git clone --recursive https://github.com/openMVG/openMVG.git
 $ cd openMVG
 $ ls
  AUTHORS BUILD  docs  logo  README  src  ...
 $ cd ..
 $ mkdir openMVG_Build
 $ cd openMVG_Build
 $ cmake -DCMAKE_BUILD_TYPE=RELEASE -G "Xcode" . ../openMVG/src/
If you want enable unit tests and examples to the build:
 $ cmake -DCMAKE_BUILD_TYPE=RELEASE -DOpenMVG_BUILD_TESTS=ON -DOpenMVG_BUILD_EXAMPLES=ON -G "Xcode" . ../openMVG/src/
 $ xcodebuild -configuration Release

 
 --------------------
 Using openCV sample
 --------------------
 Add -DUSE_OPENCV=ON to your cmake command line and set OpenCV_DIR variable to your openCV build directory
 
