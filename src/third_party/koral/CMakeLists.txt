project(koral)
cmake_minimum_required (VERSION 3.1)

if (OpenMVG_USE_CUDA)
  find_package(CUDA REQUIRED)
  find_package(SSE REQUIRED)

  include_directories(${CMAKE_CURRENT_BINARY_DIR}/config/)

  set (HEADERS CLATCH.h CUDALERP.h FeatureAngle.h KFAST.h KORAL.h Keypoint.h)

  set (SRC CLATCH.cu CUDALERP.cu)

  list(APPEND CUDA_NVCC_FLAGS "-D_MWAITXINTRIN_H_INCLUDED -lineinfo -Xptxas -v -Xcompiler -use_fast_math -O3 -gencode arch=compute_30,code=sm_30 -gencode arch=compute_60,code=sm_60 --default-stream per-thread")

  set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -std=c++11")

  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set (CMAKE_CXX_FLAGS "--std=gnu++11 ${CMAKE_CXX_FLAGS}")
  else()
    set (CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")
  endif()

  set(CUDA_PROPOGATE_HOST_FLAGS OFF)

  include_directories(${CUDA_INCLUDE_DIRS})

  cuda_add_library(koral STATIC ${SRC} ${HEADERS})
  target_link_libraries(koral ${CUDA_LIBRARIES})

  set_property(TARGET koral PROPERTY FOLDER OpenMVG/3rdParty)

  target_compile_features(koral PRIVATE cxx_constexpr)

  if(KORAL_INCLUDE_INSTALL_DIR)
    set(INCLUDE_INSTALL_DIR ${KORAL_INCLUDE_INSTALL_DIR})
  else()
    set(INCLUDE_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/include/koral")
  endif()

  install(FILES
    ${HEADERS}
    DESTINATION ${INCLUDE_INSTALL_DIR} COMPONENT Devel
  )
  install(TARGETS koral DESTINATION lib EXPORT openMVG-targets)
endif()
