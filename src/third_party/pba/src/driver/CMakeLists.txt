# Driver:

#CPU one
add_executable(driver driver.cpp)
TARGET_LINK_LIBRARIES(driver pbaCPU)

if (CUDA_FOUND)
CUDA_ADD_EXECUTABLE(driverGPU driver.cpp)
TARGET_LINK_LIBRARIES(driverGPU pbaGPU)
endif(CUDA_FOUND)