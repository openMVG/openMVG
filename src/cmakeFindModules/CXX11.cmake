if (NOT CMAKE_CXX_COMPILER_LOADED)
    message(FATAL_ERROR "CheckCXX11Features modules only works if language CXX is enabled")
endif ()

# Determines whether or not the compiler supports C++11
macro(check_for_cxx11_compiler _VAR)
  message(STATUS "Checking for C++11 compiler")
  set(${_VAR})
  if((MSVC AND NOT MSVC_VERSION VERSION_LESS 1800) OR # checking the compiler is at least Visual Studio 2013 - MSVC++ 12
     (CMAKE_COMPILER_IS_GNUCXX AND NOT ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 4.6) OR
     (CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND NOT ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 3.1))
    set(${_VAR} 1)
    message(STATUS "Checking for C++11 compiler - available")
  else()
    message(STATUS "Checking for C++11 compiler - unavailable")
  endif()
endmacro()

# Sets the appropriate flag to enable C++11 support
macro(enable_cxx11)
  if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR CMAKE_COMPILER_IS_GNUCXX)
    include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag(--std=c++11 SUPPORTS_STD_CXX11)
    check_cxx_compiler_flag(--std=c++0x SUPPORTS_STD_CXX01)
    if(SUPPORTS_STD_CXX11)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --std=c++11")
    elseif(SUPPORTS_STD_CXX01)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --std=c++0x")
    else()
      message(ERROR "Compiler does not support --std=c++11 or --std=c++0x.")
    endif()
  endif()
endmacro()
