#ifndef OPENMVG_VOCABULARY_TREE_DISTANCE_HPP
#define OPENMVG_VOCABULARY_TREE_DISTANCE_HPP

#include <stdint.h>
//#include <iostream>
#include <Eigen/Core>

namespace openMVG {
namespace voctree {

		/**
		 * \brief Meta-function returning a type that can be used to accumulate many values of T.
		 *
		 * By default, the accumulator type is the same as \c T. Specializations for the basic types
		 * are:
		 * \li \c uint8_t -> \c uint32_t
		 * \li \c uint16_t -> \c uint32_t
		 * \li \c int8_t -> \c int32_t
		 * \li \c int16_t -> \c int32_t
		 * \li \c float -> \c double
		 */
		template<typename T>
		struct Accumulator
		{
			typedef T type;
		};

		/// \cond internal

		template<> struct Accumulator<uint8_t>
		{
			typedef uint32_t type;
		};

		template<> struct Accumulator<uint16_t>
		{
			typedef uint32_t type;
		};

		template<> struct Accumulator<int8_t>
		{
			typedef int32_t type;
		};

		template<> struct Accumulator<int16_t>
		{
			typedef int32_t type;
		};
		// Have noticed loss of precision when computing sum-squared-error with many input features.

		template<> struct Accumulator<float>
		{
			typedef double type;
		};
		/// \endcond

		/**
		 * \brief Default implementation of L2 distance metric.
		 *
		 * Works with std::vector, boost::array, or more generally any container that has
		 * a \c value_type typedef, \c size() and array-indexed element access.
		 */
		template<class Feature>
		struct L2
		{
			typedef typename Feature::value_type value_type;
			typedef typename Accumulator<value_type>::type result_type;

			result_type operator( )( const Feature& a, const Feature& b ) const
			{
				result_type result = result_type( );
				for ( std::size_t i = 0; i < a.size( ); ++i )
				{
					result_type diff = a[i] - b[i];
					result += diff*diff;
				}
				return result;
			}
		};
		/// @todo Version for raw data pointers that knows the size of the feature
		/// @todo Specialization for cv::Vec. Doesn't have size() so default won't work.

		/// Specialization for Eigen::Matrix types.

		template<typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
		struct L2< Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols> >
		{
			typedef Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols> feature_type;
			typedef Scalar value_type;
			typedef typename Accumulator<Scalar>::type result_type;

			result_type operator( )( const feature_type& a, const feature_type& b ) const
			{
				//	  std::cout << a << std::endl;
				return (a - b ).squaredNorm( );
			}
		};

}
} 

#endif //OPENMVG_VOCABULARY_TREE_DISTANCE_HPP
