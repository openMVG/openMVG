set(test_target_sources)

list(APPEND test_target_sources
        QR.cpp
        Eigen.cpp
        BKLDLT.cpp
        SymEigs.cpp
        SymEigsShift.cpp
        GenEigs.cpp
        GenEigsRealShift.cpp
        GenEigsComplexShift.cpp
        SymGEigsCholesky.cpp
        SymGEigsRegInv.cpp
        SVD.cpp
        )

foreach(TEST_SOURCE ${test_target_sources})
    # Extract the filename without extension (NAME_WE) as a name for our executable
    get_filename_component(TEST_NAME ${TEST_SOURCE} NAME_WE)

    # Add an executable based on the source
    add_executable(${TEST_NAME} ${TEST_SOURCE})

    # Configure (include headers and link libraries) the test
    target_link_libraries(${TEST_NAME} PRIVATE Spectra)

    # Additional compiler flags
    target_compile_options(${TEST_NAME} PRIVATE -std=c++11 -Wall -O2 -Wno-parentheses
            )
    target_compile_options(${TEST_NAME} PRIVATE $<$<CXX_COMPILER_ID:GNU>:-Wno-misleading-indentation -Wno-int-in-bool-context>
            )

    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME}
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR})  # the working directory is the out-of-source build directory

endforeach()
