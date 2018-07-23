Multicore Bundle Adjustment

http://grail.cs.washington.edu/projects/mcba/

University of Washington at Seattle

We present the design and implementation of new inexact Newton type Bundle 
Adjustment algorithms that exploit hardware parallelism for efficiently solving 
large scale 3D scene reconstruction problems. We explore the use of  multicore 
CPU as well as multicore GPUs for this purpose. We show that overcoming the 
severe memory and bandwidth limitations of current generation GPUs not only 
leads to more space efficient algorithms, but also to surprising savings in runtime. 
Our CPU based system is up to ten times and our GPU based system is up to thirty 
times faster than the current state of the art methods, while maintaining 
comparable convergence behavior.

####################################################################
The source code is released under GUN GENERAL PUBLIC LICENSE V3. 
Please contact ccwu@cs.washington.edu for commercial licensing.
#################################################################