project(stlplus)

file(GLOB LIBSLTPLUS_HPP *.hpp )
file(GLOB LIBSLTPLUS_CPP *.cpp )

add_library(stlplus ${LIBSLTPLUS_HPP} ${LIBSLTPLUS_CPP})
