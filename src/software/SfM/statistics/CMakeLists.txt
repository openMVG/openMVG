set(OPENMVG_STATISTICS_IN_FILE "${CMAKE_CURRENT_LIST_DIR}/FeatureCountHistogram.py.in")
configure_file(${OPENMVG_STATISTICS_IN_FILE}
    "${CMAKE_CURRENT_BINARY_DIR}/FeatureCountHistogram.py" @ONLY)
