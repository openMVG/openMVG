*******************
features
*******************

This module provides generic container for features and associated descriptors.

KeypointSet 
=============

Store a collection of features and their associated descriptors: ``template<typename FeaturesT, typename DescriptorsT> class KeypointSet``.

.. code-block:: c++ 

  // The descriptor type
  typedef Descriptor<float, 128> DescriptorT;
  // The feature type
  typedef SIOPointFeature FeatureT;
  // The type of the collection to describe an image (vector):
  typedef std::vector<FeatureT> FeatsT;
  typedef std::vector<DescriptorT > DescsT;
  // A container to describe a collection of features and their descriptors:
  typedef KeypointSet<FeatsT, DescsT > KeypointSetT;


Features 
=============

Classes to store point characteristics:

 * all must inheritate from ``FeatureBase`` in order to implement serialization.

 * ``class PIOPointFeature : public PointFeature``
    * Store the position of a feature (x,y).

 * ``class SIOPointFeature : public PointFeature``
    * Store the position, orientation and scale of a feature (x,y,s,o).

Descriptors 
=============

Class to store a region description (a descriptor):

 * ``template <typename T, std::size_t N> class Descriptor``
    * Class that handle a descriptor (a data container of N values of type T).

.. code-block:: c++ 

  typedef Descriptor<float, 128> siftDescriptorData;
  typedef Descriptor<float, 64> surfDescriptorData;
 

