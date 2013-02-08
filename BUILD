=====================================
OpenMVG (open Multiple View Geometry)
=====================================

-----------------
Build instruction
-----------------

Required tools:
* Cmake 

On linux:
* => zlib, png, jpeg
* sudo apt-get install zlib1-dev libpng-dev libjpeg8-dev

-----------------
Linux compilation
-----------------

 $ ls
  AUTHORS BUILD  docs  logo  README  src  ...
 $ cd ..
 $ mkdir openMVG_Build
 $ cd openMVG_Build
 $ cmake -DCMAKE_BUILD_TYPE=RELEASE . ../openmvg/src/

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

 $ ls
  AUTHORS BUILD  docs  logo  README  src  ...
 $ cd ..
 $ mkdir openMVG_Build
 $ cd openMVG_Build
 $ cmake -DCMAKE_BUILD_TYPE=RELEASE -G "Xcode" . ../openmvg/src/
 $ xcodebuild -configuration Release

