macro(GET_OS_INFO)
    string(REGEX MATCH "Linux" OS_IS_LINUX ${CMAKE_SYSTEM_NAME})
    set(FLANN_LIB_INSTALL_DIR "lib")
    if(NOT FLANN_INCLUDE_INSTALL_DIR)
      set(FLANN_INCLUDE_INSTALL_DIR
        "include/${PROJECT_NAME_LOWER}-${FLANN_VERSION}")
    endif()    
endmacro(GET_OS_INFO)


macro(DISSECT_VERSION)
    # Find version components
    string(REGEX REPLACE "^([0-9]+).*" "\\1"
        FLANN_VERSION_MAJOR "${FLANN_VERSION}")
    string(REGEX REPLACE "^[0-9]+\\.([0-9]+).*" "\\1"
        FLANN_VERSION_MINOR "${FLANN_VERSION}")
    string(REGEX REPLACE "^[0-9]+\\.[0-9]+\\.([0-9]+)" "\\1"
        FLANN_VERSION_PATCH ${FLANN_VERSION})
    string(REGEX REPLACE "^[0-9]+\\.[0-9]+\\.[0-9]+(.*)" "\\1"
        FLANN_VERSION_CANDIDATE ${FLANN_VERSION})
    set(FLANN_SOVERSION "${FLANN_VERSION_MAJOR}.${FLANN_VERSION_MINOR}")
endmacro(DISSECT_VERSION)


# workaround a FindHDF5 bug
macro(find_hdf5)
    find_package(HDF5)

    set( HDF5_IS_PARALLEL FALSE )
    foreach( _dir ${HDF5_INCLUDE_DIRS} )
        if( EXISTS "${_dir}/H5pubconf.h" )
            file( STRINGS "${_dir}/H5pubconf.h" 
                HDF5_HAVE_PARALLEL_DEFINE
                REGEX "HAVE_PARALLEL 1" )
            if( HDF5_HAVE_PARALLEL_DEFINE )
                set( HDF5_IS_PARALLEL TRUE )
            endif()
        endif()
    endforeach()
    set( HDF5_IS_PARALLEL ${HDF5_IS_PARALLEL} CACHE BOOL
        "HDF5 library compiled with parallel IO support" )
    mark_as_advanced( HDF5_IS_PARALLEL )
endmacro(find_hdf5)


macro(flann_add_gtest exe)
    # add build target
    add_executable(${exe} EXCLUDE_FROM_ALL ${ARGN})
    target_link_libraries(${exe} ${GTEST_LIBRARIES})
    # add dependency to 'tests' target
    add_dependencies(flann_gtests ${exe})

    # add target for running test
    string(REPLACE "/" "_" _testname ${exe})
    add_custom_target(test_${_testname}
                    COMMAND ${exe}
                    ARGS --gtest_print_time
                    DEPENDS ${exe}
                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test
                    VERBATIM
                    COMMENT "Runnint gtest test(s) ${exe}")
    # add dependency to 'test' target
    add_dependencies(flann_gtest test_${_testname})
endmacro(flann_add_gtest)

macro(flann_add_cuda_gtest exe)
    # add build target
    cuda_add_executable(${exe} EXCLUDE_FROM_ALL ${ARGN})
    target_link_libraries(${exe} ${GTEST_LIBRARIES})
    # add dependency to 'tests' target
    add_dependencies(tests ${exe})

    # add target for running test
    string(REPLACE "/" "_" _testname ${exe})
    add_custom_target(test_${_testname}
                    COMMAND ${exe}
                    ARGS --gtest_print_time
                    DEPENDS ${exe}
                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test
                    VERBATIM
                    COMMENT "Runnint gtest test(s) ${exe}")
    # add dependency to 'test' target
    add_dependencies(test test_${_testname})
endmacro(flann_add_cuda_gtest)

macro(flann_add_pyunit file)
    # find test file
    set(_file_name _file_name-NOTFOUND)
    find_file(_file_name ${file} ${CMAKE_CURRENT_SOURCE_DIR})
    if(NOT _file_name)
        message(FATAL_ERROR "Can't find pyunit file \"${file}\"")
    endif(NOT _file_name)

    # add target for running test
    string(REPLACE "/" "_" _testname ${file})
    add_custom_target(pyunit_${_testname}
                    COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_SOURCE_DIR}/bin/run_test.py ${_file_name}
                    DEPENDS ${_file_name}
                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test
                    VERBATIM
                    COMMENT "Running pyunit test(s) ${file}" )
    # add dependency to 'test' target
    add_dependencies(pyunit_${_testname} flann)
    add_dependencies(test pyunit_${_testname})
endmacro(flann_add_pyunit)



macro(flann_download_test_data _name _md5)
    string(REPLACE "/" "_" _dataset_name dataset_${_name})
    
    add_custom_target(${_dataset_name}
        COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_SOURCE_DIR}/bin/download_checkmd5.py http://people.cs.ubc.ca/~mariusm/uploads/FLANN/datasets/${_name} ${TEST_OUTPUT_PATH}/${_name} ${_md5}
        VERBATIM)

    # Also make sure that downloads are done before we run any tests
    add_dependencies(tests ${_dataset_name})

endmacro(flann_download_test_data)
