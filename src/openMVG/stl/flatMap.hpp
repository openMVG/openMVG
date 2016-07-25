#pragma once

#ifdef HAVE_BOOST
#include <boost/container/flat_map.hpp>
#include <cereal/types/map.hpp>

// Add serialization for boost flat_map
namespace cereal
{
  //! Saving for std::map
  template <class Archive, class K, class T, class C, class A> inline
  void CEREAL_SAVE_FUNCTION_NAME( Archive & ar, boost::container::flat_map<K, T, C, A> const & map )
  {
    map_detail::save( ar, map );
  }

  //! Loading for std::map
  template <class Archive, class K, class T, class C, class A> inline
  void CEREAL_LOAD_FUNCTION_NAME( Archive & ar, boost::container::flat_map<K, T, C, A> & map )
  {
    map_detail::load( ar, map );
  }
} // namespace cereal
#endif

namespace stl
{
#ifdef HAVE_BOOST
  template <class Key
         ,class T
         ,class Compare = std::less<Key>
         ,class Allocator = std::allocator<std::pair<Key, T> > >
  using flat_map = boost::container::flat_map<Key, T, Compare, Allocator>;
#else
  // Fallback to non-flat implementation
  template <class Key
         ,class T
         ,class Compare = std::less<Key>
         ,class Allocator = std::allocator<std::pair<Key, T> > >
  using flat_map = std::map<Key, T, Compare, Allocator>;
#endif

}
