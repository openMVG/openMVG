#ifndef OPENMVG_STL_BITMASK_H
#define OPENMVG_STL_BITMASK_H

// Taken from boost
#define OPENMVG_BITMASK(Bitmask)                                            \
                                                                          \
  inline Bitmask operator| (Bitmask x , Bitmask y )                       \
  { return static_cast<Bitmask>( static_cast<std::int32_t>(x)     \
      | static_cast<std::int32_t>(y)); }                          \
                                                                          \
  inline Bitmask operator& (Bitmask x , Bitmask y )                       \
  { return static_cast<Bitmask>( static_cast<std::int32_t>(x)     \
      & static_cast<std::int32_t>(y)); }                          \
                                                                          \
  inline Bitmask operator^ (Bitmask x , Bitmask y )                       \
  { return static_cast<Bitmask>( static_cast<std::int32_t>(x)     \
      ^ static_cast<std::int32_t>(y)); }                          \
                                                                          \
  inline Bitmask operator~ (Bitmask x )                                   \
  { return static_cast<Bitmask>(~static_cast<std::int32_t>(x)); } \
                                                                          \
  inline Bitmask & operator&=(Bitmask & x , Bitmask y)                    \
  { x = x & y ; return x ; }                                              \
                                                                          \
  inline Bitmask & operator|=(Bitmask & x , Bitmask y)                    \
  { x = x | y ; return x ; }                                              \
                                                                          \
  inline Bitmask & operator^=(Bitmask & x , Bitmask y)                    \
  { x = x ^ y ; return x ; }                                              


#endif  // OPENMVG_STL_HASH_H
