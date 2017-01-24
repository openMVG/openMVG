*******************
features
*******************

This module provides generic container for features and associated descriptors.

Features
=============

Provide basic structure and IO to store Point based features.

Classes to store point characteristics:

 * :class:`PointFeature`
    * Store the position of a feature (x,y).

 * :class:`SIOPointFeature`
    * Store the position, orientation and scale of a feature (x,y,s,o).

Descriptors
=============

Provide basic structure and IO for descriptor data.

 * :class:`template <typename T, std::size_t N> class Descriptor`.
    * Store N value(s) of type T as contiguous memory.

.. code-block:: c++

  // SIFT like descriptor
  using siftDescriptorData Descriptor<float, 128>;

  // SURF like descriptor
  using surfDescriptorData = Descriptor<float, 64>;

  // Binary descriptor (128 bits)
  using binaryDescriptor_bitset = Descriptor<std::bitset<128>,1> binaryDescriptor_bitset;
  // or using unsigned chars
  using binaryDescriptor_uchar = Descriptor<unsigned char, 128/sizeof(unsigned char)>;


KeypointSet
=============

Store a collection of features and their associated descriptors: :class:`template<typename FeaturesT, typename DescriptorsT> class KeypointSet`. Basic IO is provided.

.. code-block:: c++

  // Define SIFT Keypoints:

  // Define the SIFT descriptor [128 floating point value]
  using Descriptor<float, 128> DescriptorT;

  // Use SIFT compatible features (scale, orientation and position)
  using FeatureT = SIOPointFeature;

  // Describe what a collection of local feature is for a given image:
  using FeatsT = std::vector<FeatureT>;
  using DescsT = std::vector<DescriptorT>;

  // Link features and their descriptors as a collection:
  using KeypointSetT = KeypointSet<FeatsT, DescsT>;
