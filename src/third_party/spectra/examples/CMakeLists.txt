set(example_target_sources)

list(APPEND example_target_sources
        DavidsonSymEigs_example.cpp
        )

foreach(EXAMPLE_SOURCE ${example_target_sources})
    # Extract the filename without extension (NAME_WE) as a name for our executable
    get_filename_component(EXAMPLE_NAME ${EXAMPLE_SOURCE} NAME_WE)

    # Add an executable based on the source
    add_executable(${EXAMPLE_NAME} ${EXAMPLE_SOURCE})

    # Configure (include headers and link libraries) the example
    target_link_libraries(${EXAMPLE_NAME} PRIVATE Spectra)

endforeach()
