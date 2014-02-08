*******************
multiview
*******************

The multiview module consists of a collection of:

- solvers for 2 to n-view geometry constraints that arise in multiple view geometry,
- a generic framework "Kernel" that can embed these solvers for robust estimation.

First accessible solvers are listed and explained and the "Kernel" concept is documented.

Two-view 2d-2d geometric estimation
=====================================

openMVG provides solver for the following problems:

* affine,
*	homographic,
* fundamental,

  * 7 to n pt,
  * 8 to n pt (Direct Linear Transform) [HZ]_.
  
* essential,

  * 8 to n pt (Direct Linear Transform) [HZ]_,
  * 5pt + intrinsic [Stewenius]_, [Nister]_.
    
One-view 2d-3d geometric estimation
=====================================

openMVG provides solver for the following problems:

* pose estimation/ camera resection,

  * 6pt Direct Linear Transform [HZ]_,
  * 4pt with intrinsic EPnP [EPnP]_,  
  * 3pt with intrinsic P3P [Kneip]_.


N-View geometry estimation
============================

* Triangulation
  * 2 to n view (Direct Linear Transform),
  * 2 to n view (Iterated least square).
 
Homography matrix:
---------------------

The homography matrix maps the relation between two projections of a plane: :ref:`Figure <homographyFigure>`.

:math:`H` is a (3 x 3) matrix that links coordinates in left and right images with the following relation.

.. math:: x'_i = Hx_i
   :label: homographyformula

.. _homographyFigure:
.. figure:: homographyMatrix.png
   :align: center

   The homography matrix and the point to point constraint.

OpenMVG implementation follows the DLT (Direct Linear Transform) explained in [HZ]_ book:
:math:`H` can be estimated from 4 to :math:`n` corresponding points.
   
Fundamental matrix:
---------------------

The fundamental matrix is a relation between two images viewing the same scene where those point's
projections are visible in the two images. Given a point correspondence between two views :math:`(x_i,x'_i)`:

We obtain the following relation:

.. math:: x'^T_i Fx_i = 0

:math:`F` is the (3 x 3) Fundamental matrix, it puts in relation a point x to a line where belong the projection of the 3D X point.
:math:`l'_i = Fx_i` designs the epipolar line on which the point :math:`x'_i` could be.
The relation :math:`x'^T_i Fx_i = 0` exists for all corresponding point belonging to a stereo pair.

.. figure:: fundamentalMatrix.png
   :align: center

   The fundamental matrix and the point to line constraint.

The fundamental matrix is sometime called bifocal-tensor, it is a 3 x 3 matrix of rank 2
with 7 degree of freedom. 8 ou 7 correspondences are sufficient to compute the :math:`F` matrix.
Implementation follows the DLT (Direct Linear Transform) explained in [HZ]_ book.

Essential matrix
---------------------

Adding intrinsic parameters to the fundamental matrix gives a metric "object" that provides the following relation
:math:`E = K'^T FK` , this is the Essential relation explained by Longuet-Higgins in 1981 [Longuet]_.
This essential matrix links the relative position of the camera to the fundamental matrix relation.

.. math:: E = R[t]x = K'^T FK

.. figure:: essentialMatrix.png
   :align: center

   The essential matrix geometric relation.

Camera resection/pose estimation
===================================

Given a serie of 3D-2D image plane correspondences it is possible to compute a camera pose estimation.
It consists in estimating the camera parameters of the right camera that minimizes the residual error of the 3D points re-projections, it's an optimization problem that trying to solve P parameters in order to minimize 

.. math:: 
	min \sum^n_{i=1} x_i - P(X_i).


.. figure:: ResectionPose.png
   :align: center
	
   Resection/Pose estimation from 3D-2D correspondences.

.. figure:: ResidualError.png
   :align: center
	
   Residual error.

Kernel concept
---------------------

A kernel is an association :

* data points (the set that is used for a robust estimation problem),
* a model solver/estimator
* a metric to measure data fitting to a putative model.

.. figure:: kernelConcept.png
   :align: center

   The Kernel concept (the two view case).

Solver:
	MINIMUM_SAMPLES: The minimal number of point required for the model estimation,
	MAX_MODELS: The number of models that the minimal solver could return,
	A Solve function that estimates a model from MINIMUM_SAMPLES to n vector data.
ErrorModel: 
	An metric function that return the error of a sample data to the provided model.
Kernel: 
	Embed data (putative), the model estimator and the error model.
	This kernel is core brick used in the openMVG robust estimator framework.

.. [Kneip] A Novel Parametrization of the P3P-Problem for a Direct Computation of Absolute Camera Position and Orientation.
    Kneip, L.; Scaramuzza, D. ; Siegwart, R.
    CVPR 2011
    
.. [EPnP] EPnP: An Accurate O(n) Solution to the PnP Problem.
    V. Lepetit and F. Moreno-Noguer and P. Fua, IJCV 2009. vol. 81, no. 2.
    
.. [HZ] Multiple view geometry in computer vision.
    Hartley, Richard, and Andrew Zisserman. 
    Vol. 2. Cambridge, 2000.
    
.. [Stewenius] Recent Developments on Direct Relative Orientation.
      H. Stewenius, C. Engels and D. Nister.
      ISPRS 2006

.. [Nister] An Efficient Solution to the Five-Point Relative Pose.
    D. Nister
    PAMI 2004
    
.. [Longuet] A computer algorithm for reconstructing a scene from two projections.
    Longuet-Higgins, H. C.
    Readings in Computer Vision: Issues, Problems, Principles, and Paradigms, MA Fischler and O. Firschein, eds (1987): 61-62.
