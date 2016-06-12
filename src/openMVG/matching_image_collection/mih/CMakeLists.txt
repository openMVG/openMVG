cmake_minimum_required(VERSION 2.6)
project (mih_cpu)

find_library (HDF5_LIB hdf5)

set (CMAKE_CXX_FLAGS "-DSINGLE_CORE -O2 -fopenmp -DNDEBUG")
# -Wall -march=native -mtune=native")

include_directories("./include")

add_library (mih_cpu src/array32.cpp src/mihasher.cpp src/sparse_hashtable.cpp src/bucket_group.cpp src/reorder.cpp
					 include/array32.h include/bitops.h include/io.h include/memusage.h include/myhdf5.h include/result.h
					 include/sparse_hashtable.h include/bitarray.h include/bucket_group.h include/linscan.h 
					 include/mihasher.h include/reorder.h include/resulth5.h include/types.h)
target_link_libraries (mih_cpu ${HDF5_LIB})

install(TARGETS mih_cpu DESTINATION lib EXPORT openMVG-targets)
