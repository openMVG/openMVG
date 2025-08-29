/* -*- mode: C++; indent-tabs-mode: nil; -*-
 */

#ifndef STL_ITERATORS_H_
#define STL_ITERATORS_H_

#include <lemon/core.h>

namespace lemon {

  /// \brief Template to make STL iterators from Lemon iterators.
  ///
  /// This template makes an STL iterator from a Lemon iterator
  /// by adding the missing features.
  /// It inherits from \c std::iterator to make \c iterator_concept work
  /// (so that STL algorithms work).
  /// \c T should be the lemon iterator to be decorated.
  template<class T>
  struct LemonItWrapper
      : public T, public std::iterator<std::input_iterator_tag, T> {

    LemonItWrapper(const T &x) : T(x) {}

    //Lemon iterators don't have operator*, (because they rather
    //inherit from their "value_type"),
    //so we add one that just returns the object.
    const T& operator*() const {
      return static_cast<const T&>(*this);
    }

    //I can't think of any use case for this with Lemon iterators,
    //but maybe it should be included for completeness.
    const T* operator->() {
      return static_cast<const T*>(this);
    }

    //Lemon iterators don't have postincrement.
    void operator++(int) {
      T::operator++();
    }

    using T::operator++;

  };


  /// \brief A generic wrapper for Lemon iterators for range-based for loops.
  ///
  /// This template can be used to create a class
  /// that has begin() and end() from a Lemon iterator
  /// (with a 1-parameter constructor)
  /// to make range-based for loops and STL algorithms work.
  ///
  /// \c LIT is the Lemon iterator that will be wrapped
  /// \c P is the type of the parameter of the constructor of \c LIT.
  template<class LIT, class P>
  class LemonRangeWrapper1 {
    typedef LemonItWrapper<LIT> It;
    It _begin;
    
  public:
    LemonRangeWrapper1(const P &p) : _begin(LIT(p)) {}
    It begin() const {
      return _begin;
    }
    It end() const {
      return It(lemon::INVALID);
    }
  };


  /// \brief A generic wrapper for Lemon iterators for range-based for loops.
  ///
  /// This template can be used to create a class
  /// that has begin() and end() from a Lemon iterator
  /// (with a 2-parameter constructor)
  /// to make range-based for loops and STL algorithms work.
  ///
  /// \c LIT is the Lemon iterator that will be wrapped
  /// \c P1 and \c P2 are the types of the parameters
  /// of the constructor of \c LIT.
  template<class LIT, class P1, class P2>
  class LemonRangeWrapper2 {
    typedef LemonItWrapper<LIT> It; 
    It _begin;
 public:
    LemonRangeWrapper2(const P1 &p1, const P2 &p2) : _begin(LIT(p1, p2)) {}
    It begin() const {
      return _begin;
    }
    It end() const {
      return It(lemon::INVALID);
    }
  };


}

#endif /* STL_ITERATORS_H_ */
