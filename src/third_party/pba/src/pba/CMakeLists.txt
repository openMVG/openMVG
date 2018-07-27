
SET(PBA_CPU_SRCS
    ConfigBA.h
    ConfigBA.cpp
    DataInterface.h
    interface.h
    pba.h
    pba.cpp
    SparseBundleCPU.h
    SparseBundleCPU.cpp
    util.h)

add_library(pbaCPU ${PBA_CPU_SRCS})
SET_TARGET_PROPERTIES(pbaCPU PROPERTIES COMPILE_FLAGS "-DPBA_NO_GPU")

if (CUDA_FOUND) 
  SET(PBA_GPU_SRCS
      ConfigBA.h
      ConfigBA.cpp
      DataInterface.h
      interface.h
      pba.h
      pba.cpp
      SparseBundleCU.h
      SparseBundleCU.cpp
      SparseBundleCPU.h
      SparseBundleCPU.cpp
      CuTexImage.h
      CuTexImage.cpp
      ProgramCU.h
      ProgramCU.cu
      util.h)

  CUDA_ADD_LIBRARY(pbaGPU ${PBA_GPU_SRCS})
endif(CUDA_FOUND) 

install(TARGETS pbaCPU DESTINATION lib EXPORT openMVG-targets)