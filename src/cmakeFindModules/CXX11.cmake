if (NOT CMAKE_CXX_COMPILER_LOADED)
    message(FATAL_ERROR "CheckCXX11Features modules only works if language CXX is enabled")
endif ()

# Determines whether or not the compiler supports C++17
# (project default standard is C++17)
macro(check_for_cxx17_compiler _VAR)
  message(STATUS "Checking for C++17 compiler")
  set(${_VAR})
  if((MSVC AND NOT MSVC_VERSION VERSION_LESS 1914) OR # checking the compiler is at least Visual Studio 2017 - MSVC++ 14.1
     (CMAKE_COMPILER_IS_GNUCXX AND NOT ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 7.0) OR
     (CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND NOT ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 5.0))
    set(${_VAR} 1)
    message(STATUS "Checking for C++17 compiler - available")
  else()
    message(STATUS "Checking for C++17 compiler - unavailable")
  endif()
endmacro()

if(NOT (${CMAKE_VERSION} VERSION_LESS "3.8.0"))
  # For CMake 3.8 and above, we can use meta features directly provided by CMake itself
  set(CXX11_FEATURES cxx_std_11)
  set(CXX14_FEATURES cxx_std_14)
  set(CXX17_FEATURES cxx_std_17)
  return()
endif()

set(CXX11_FEATURES
  cxx_auto_type
  cxx_constexpr
  cxx_lambdas
  cxx_nullptr
  cxx_override
  cxx_range_for
  cxx_strong_enums
)
