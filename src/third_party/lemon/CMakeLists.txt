CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

SET(PROJECT_NAME "LEMON")
PROJECT(${PROJECT_NAME})

INCLUDE(FindPythonInterp)
INCLUDE(FindWget)

IF(EXISTS ${PROJECT_SOURCE_DIR}/cmake/version.cmake)
  INCLUDE(${PROJECT_SOURCE_DIR}/cmake/version.cmake)
ELSEIF(DEFINED ENV{LEMON_VERSION})
  SET(LEMON_VERSION $ENV{LEMON_VERSION} CACHE STRING "LEMON version string.")
ELSE()
  EXECUTE_PROCESS(
    COMMAND ${PYTHON_EXECUTABLE} ./scripts/chg-len.py
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE HG_REVISION_PATH
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  EXECUTE_PROCESS(
    COMMAND hg id -i
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE HG_REVISION
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  IF(HG_REVISION STREQUAL "")
    SET(HG_REVISION_ID "hg-tip")
  ELSE()
    IF(HG_REVISION_PATH STREQUAL "")
      SET(HG_REVISION_ID ${HG_REVISION})
    ELSE()
      SET(HG_REVISION_ID ${HG_REVISION_PATH}.${HG_REVISION})
    ENDIF()
  ENDIF()
  SET(LEMON_VERSION ${HG_REVISION_ID} CACHE STRING "LEMON version string.")
ENDIF()

SET(PROJECT_VERSION ${LEMON_VERSION})

SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

FIND_PACKAGE(Doxygen)
FIND_PACKAGE(Ghostscript)
FIND_PACKAGE(GLPK 4.33)
FIND_PACKAGE(CPLEX)
FIND_PACKAGE(COIN)

IF(DEFINED ENV{LEMON_CXX_WARNING})
  SET(CXX_WARNING $ENV{LEMON_CXX_WARNING})
ELSE()
  IF(CMAKE_COMPILER_IS_GNUCXX)
    SET(CXX_WARNING "-Wall -W -Wunused -Wformat=2 -Wctor-dtor-privacy -Wnon-virtual-dtor -Wno-char-subscripts -Wwrite-strings -Wno-char-subscripts -Wreturn-type -Wcast-qual -Wcast-align -Wsign-promo -Woverloaded-virtual -fno-strict-aliasing -Wold-style-cast -Wno-unknown-pragmas")
    SET(CMAKE_CXX_FLAGS_DEBUG CACHE STRING "-ggdb")
    SET(CMAKE_C_FLAGS_DEBUG CACHE STRING "-ggdb")
  ELSEIF(MSVC)
    # This part is unnecessary 'casue the same is set by the lemon/core.h.
    # Still keep it as an example.
    SET(CXX_WARNING "/wd4250 /wd4355 /wd4503 /wd4800 /wd4996")
    # Suppressed warnings:
    # C4250: 'class1' : inherits 'class2::member' via dominance
    # C4355: 'this' : used in base member initializer list
    # C4503: 'function' : decorated name length exceeded, name was truncated
    # C4800: 'type' : forcing value to bool 'true' or 'false'
    #        (performance warning)
    # C4996: 'function': was declared deprecated
  ELSE()
    SET(CXX_WARNING "-Wall -W")
  ENDIF()
ENDIF()
SET(LEMON_CXX_WARNING_FLAGS ${CXX_WARNING} CACHE STRING "LEMON warning flags.")

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LEMON_CXX_WARNING_FLAGS}")

SET( CMAKE_CXX_FLAGS_MAINTAINER "-Werror -ggdb" CACHE STRING
    "Flags used by the C++ compiler during maintainer builds."
    FORCE )
SET( CMAKE_C_FLAGS_MAINTAINER "-Werror" CACHE STRING
    "Flags used by the C compiler during maintainer builds."
    FORCE )
SET( CMAKE_EXE_LINKER_FLAGS_MAINTAINER
    "-Wl,--warn-unresolved-symbols,--warn-once" CACHE STRING
    "Flags used for linking binaries during maintainer builds."
    FORCE )
SET( CMAKE_SHARED_LINKER_FLAGS_MAINTAINER
    "-Wl,--warn-unresolved-symbols,--warn-once" CACHE STRING
    "Flags used by the shared libraries linker during maintainer builds."
    FORCE )
MARK_AS_ADVANCED(
    CMAKE_CXX_FLAGS_MAINTAINER
    CMAKE_C_FLAGS_MAINTAINER
    CMAKE_EXE_LINKER_FLAGS_MAINTAINER
    CMAKE_SHARED_LINKER_FLAGS_MAINTAINER )

IF(CMAKE_CONFIGURATION_TYPES)
  LIST(APPEND CMAKE_CONFIGURATION_TYPES Maintainer)
  LIST(REMOVE_DUPLICATES CMAKE_CONFIGURATION_TYPES)
  SET(CMAKE_CONFIGURATION_TYPES "${CMAKE_CONFIGURATION_TYPES}" CACHE STRING
      "Add the configurations that we need"
      FORCE)
 endif()

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE "Release")
ENDIF()

SET( CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING
    "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel Maintainer."
    FORCE )


INCLUDE(CheckTypeSize)
CHECK_TYPE_SIZE("long long" LONG_LONG)
SET(LEMON_HAVE_LONG_LONG ${HAVE_LONG_LONG})

INCLUDE(FindPythonInterp)

ENABLE_TESTING()

IF(${CMAKE_BUILD_TYPE} STREQUAL "Maintainer")
  ADD_CUSTOM_TARGET(check ALL COMMAND ${CMAKE_CTEST_COMMAND})
ELSE()
  ADD_CUSTOM_TARGET(check COMMAND ${CMAKE_CTEST_COMMAND})
ENDIF()

ADD_SUBDIRECTORY(lemon)
IF(${CMAKE_SOURCE_DIR} STREQUAL ${PROJECT_SOURCE_DIR})
  ADD_SUBDIRECTORY(demo)
  ADD_SUBDIRECTORY(tools)
  ADD_SUBDIRECTORY(doc)
  ADD_SUBDIRECTORY(test)
ENDIF()

CONFIGURE_FILE(
  ${PROJECT_SOURCE_DIR}/cmake/LEMONConfig.cmake.in
  ${PROJECT_BINARY_DIR}/cmake/LEMONConfig.cmake
  @ONLY
)
IF(UNIX)
  INSTALL(
    FILES ${PROJECT_BINARY_DIR}/cmake/LEMONConfig.cmake
    DESTINATION share/lemon/cmake
  )
ELSEIF(WIN32)
  INSTALL(
    FILES ${PROJECT_BINARY_DIR}/cmake/LEMONConfig.cmake
    DESTINATION cmake
  )
ENDIF()

IF(${CMAKE_SOURCE_DIR} STREQUAL ${PROJECT_SOURCE_DIR})
  SET(CPACK_PACKAGE_NAME ${PROJECT_NAME})
  SET(CPACK_PACKAGE_VENDOR "EGRES")
  SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY
    "LEMON - Library for Efficient Modeling and Optimization in Networks")
  SET(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE")

  SET(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})

  SET(CPACK_PACKAGE_INSTALL_DIRECTORY
    "${PROJECT_NAME} ${PROJECT_VERSION}")
  SET(CPACK_PACKAGE_INSTALL_REGISTRY_KEY
    "${PROJECT_NAME} ${PROJECT_VERSION}")

  SET(CPACK_COMPONENTS_ALL headers library html_documentation bin)

  SET(CPACK_COMPONENT_HEADERS_DISPLAY_NAME "C++ headers")
  SET(CPACK_COMPONENT_LIBRARY_DISPLAY_NAME "Dynamic-link library")
  SET(CPACK_COMPONENT_BIN_DISPLAY_NAME "Command line utilities")
  SET(CPACK_COMPONENT_HTML_DOCUMENTATION_DISPLAY_NAME "HTML documentation")

  SET(CPACK_COMPONENT_HEADERS_DESCRIPTION
    "C++ header files")
  SET(CPACK_COMPONENT_LIBRARY_DESCRIPTION
    "DLL and import library")
  SET(CPACK_COMPONENT_BIN_DESCRIPTION
    "Command line utilities")
  SET(CPACK_COMPONENT_HTML_DOCUMENTATION_DESCRIPTION
    "Doxygen generated documentation")

  SET(CPACK_COMPONENT_HEADERS_DEPENDS library)

  SET(CPACK_COMPONENT_HEADERS_GROUP "Development")
  SET(CPACK_COMPONENT_LIBRARY_GROUP "Development")
  SET(CPACK_COMPONENT_HTML_DOCUMENTATION_GROUP "Documentation")

  SET(CPACK_COMPONENT_GROUP_DEVELOPMENT_DESCRIPTION
    "Components needed to develop software using LEMON")
  SET(CPACK_COMPONENT_GROUP_DOCUMENTATION_DESCRIPTION
    "Documentation of LEMON")

  SET(CPACK_ALL_INSTALL_TYPES Full Developer)

  SET(CPACK_COMPONENT_HEADERS_INSTALL_TYPES Developer Full)
  SET(CPACK_COMPONENT_LIBRARY_INSTALL_TYPES Developer Full)
  SET(CPACK_COMPONENT_HTML_DOCUMENTATION_INSTALL_TYPES Full)

  SET(CPACK_GENERATOR "NSIS")
  SET(CPACK_NSIS_MUI_ICON "${PROJECT_SOURCE_DIR}/cmake/nsis/lemon.ico")
  SET(CPACK_NSIS_MUI_UNIICON "${PROJECT_SOURCE_DIR}/cmake/nsis/uninstall.ico")
  #SET(CPACK_PACKAGE_ICON "${PROJECT_SOURCE_DIR}/cmake/nsis\\\\installer.bmp")
  SET(CPACK_NSIS_INSTALLED_ICON_NAME "bin\\\\lemon.ico")
  SET(CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY} ${PROJECT_NAME}")
  SET(CPACK_NSIS_HELP_LINK "http:\\\\\\\\lemon.cs.elte.hu")
  SET(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\lemon.cs.elte.hu")
  SET(CPACK_NSIS_CONTACT "lemon-user@lemon.cs.elte.hu")
  SET(CPACK_NSIS_CREATE_ICONS_EXTRA "
    CreateShortCut \\\"$SMPROGRAMS\\\\$STARTMENU_FOLDER\\\\Documentation.lnk\\\" \\\"$INSTDIR\\\\share\\\\doc\\\\index.html\\\"
    ")
  SET(CPACK_NSIS_DELETE_ICONS_EXTRA "
    !insertmacro MUI_STARTMENU_GETFOLDER Application $MUI_TEMP
    Delete \\\"$SMPROGRAMS\\\\$MUI_TEMP\\\\Documentation.lnk\\\"
    ")

  INCLUDE(CPack)
ENDIF()
