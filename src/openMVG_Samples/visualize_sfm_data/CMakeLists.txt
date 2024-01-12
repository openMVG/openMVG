#rerun_sdk
include(FetchContent)
FetchContent_Declare(rerun_sdk URL
    https://github.com/rerun-io/rerun/releases/latest/download/rerun_cpp_sdk.zip)
FetchContent_MakeAvailable(rerun_sdk)

add_executable(openMVG_sample_visualize_sfm_data visualize_sfm_data.cpp)
target_link_libraries(openMVG_sample_visualize_sfm_data PRIVATE rerun_sdk openMVG_sfm openMVG_image ${STLPLUS_LIBRARY})
target_compile_features(openMVG_sample_visualize_sfm_data PRIVATE cxx_std_17)
target_compile_definitions(openMVG_sample_visualize_sfm_data PRIVATE -DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
