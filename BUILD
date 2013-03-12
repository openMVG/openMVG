=====================================
OpenMVG (open Multiple View Geometry)
=====================================

-----------------
Build instruction
-----------------

Required tools:
* Cmake 
* Git

Getting the sources:
$ git clone --recursive https://github.com/openMVG/openMVG.git

Depending of your platform :

-----------------
Linux compilation
-----------------

Setup the required external library.
* sudo apt-get install libpng-dev libjpeg-dev x11proto-xf86vidmode-dev libxrandr-dev

 $ git clone --recursive https://github.com/openMVG/openMVG.git
 $ cd openMVG
 $ ls
  AUTHORS BUILD  docs  logo  README  src  ...
 $ cd ..
 $ mkdir openMVG_Build
 $ cd openMVG_Build
 $ cmake -DCMAKE_BUILD_TYPE=RELEASE . ../openMVG/src/

If you want have an IDE openable project with codeblocks:
 $ cmake -G "CodeBlocks - Unix Makefiles" -DCMAKE_BUILD_TYPE=RELEASE . ../openmvg/src/

Compile the project
 $ make

For a multi-core compilation (Replace NBcore with the number of threads)
 $ make -j NBcore

Launch test
 $ make test

Have fun with the samples
 $ cd opengMVG_Samples

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
 $ cmake -DCMAKE_BUILD_TYPE=RELEASE -G "Xcode" . ../openmvg/src/
 $ xcodebuild -configuration Release

