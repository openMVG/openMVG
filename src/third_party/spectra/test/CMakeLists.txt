set(test_target_sources)


add_library (tests-main tests-main.cpp)
list(APPEND test_target_sources
        BKLDLT.cpp
        DavidsonSymEigs.cpp
        DenseGenMatProd.cpp
        DenseSymMatProd.cpp
        Eigen.cpp
        GenEigs.cpp
        GenEigsRealShift.cpp
        GenEigsComplexShift.cpp
        Orthogonalization.cpp
        JDSymEigsBase.cpp
        JDSymEigsDPRConstructor.cpp
        QR.cpp
        RitzPairs.cpp
        Schur.cpp
        SearchSpace.cpp
        SparseGenMatProd.cpp
        SparseSymMatProd.cpp
        SVD.cpp
        SymEigs.cpp
        SymEigsShift.cpp
        SymGEigsCholesky.cpp
        SymGEigsRegInv.cpp
        SymGEigsShift.cpp
        )

foreach(TEST_SOURCE ${test_target_sources})
    # Extract the filename without extension (NAME_WE) as a name for our executable
    get_filename_component(TEST_NAME ${TEST_SOURCE} NAME_WE)

    # Add an executable based on the source
    add_executable(${TEST_NAME} ${TEST_SOURCE})

    # Configure (include headers and link libraries) the test
    target_link_libraries(${TEST_NAME} PRIVATE Spectra tests-main)

    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME}
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR})  # the working directory is the out-of-source build directory
endforeach()

