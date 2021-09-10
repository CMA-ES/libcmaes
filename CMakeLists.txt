cmake_minimum_required (VERSION 3.1)

project (libcmaes VERSION 0.10 LANGUAGES C CXX)

if (NOT DEFINED CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release CACHE STRING "Build type")
endif ()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

list (APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

set (LIBCMAES_TOP_LEVEL NO)
if (${PROJECT_SOURCE_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
  set (LIBCMAES_TOP_LEVEL YES)
endif ()

# ---------- options ----------
option (LIBCMAES_BUILD_SHARED_LIBS "Build libcmaes as a shared library" ON)
option (LIBCMAES_BUILD_PYTHON "build python bindings" OFF)
option (LIBCMAES_BUILD_TESTS "build tests" OFF)
option (LIBCMAES_BUILD_EXAMPLES "build examples" ${LIBCMAES_TOP_LEVEL})
option (LIBCMAES_USE_OPENMP "Use OpenMP for multithreading" ON)
option (LIBCMAES_ENABLE_SURROG "support for surrogates" ON)

# Offer the user the choice of overriding the installation directories
set (INSTALL_LIB_DIR lib${LIB_SUFFIX}
     CACHE PATH "Installation directory for libraries")
set (INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
set (INSTALL_INCLUDE_DIR include
     CACHE PATH "Installation directory for header files")
set (INSTALL_DATA_DIR share/libcmaes
     CACHE PATH "Installation directory for data files")
set (INSTALL_CMAKE_DIR ${INSTALL_LIB_DIR}/cmake/libcmaes
     CACHE PATH "Installation directory for cmake config files")

# Make relative paths absolute (needed later on)
foreach (p LIB BIN INCLUDE DATA CMAKE)
  set (var INSTALL_${p}_DIR)
  set (RELATIVE_INSTALL_${p}_DIR ${INSTALL_${p}_DIR})
  if (NOT IS_ABSOLUTE "${${var}}")
    set (${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif ()
endforeach ()

set (CMAKE_INSTALL_RPATH ${INSTALL_LIB_DIR})
set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# ----------- required include files ----------
include (CheckIncludeFile)

check_include_file (dlfcn.h HAVE_DLFCN_H)
check_include_file (unistd.h HAVE_UNISTD_H)
check_include_file (string.h HAVE_STRING_H)
check_include_file (strings.h HAVE_STRINGS_H)
check_include_file (inttypes.h HAVE_INTTYPES_H)
check_include_file (memory.h HAVE_MEMORY_H)
check_include_file (stdlib.h HAVE_STDLIB_H)
check_include_file (stdint.h HAVE_STDINT_H)
check_include_file (sys/types.h HAVE_SYS_TYPES_H)
check_include_file (sys/stat.h HAVE_SYS_STAT_H)

# ---------- dependencies ----------
find_package (Eigen3 3.0.0 REQUIRED)

if (LIBCMAES_USE_OPENMP)
  find_package (OpenMP)
  if (NOT TARGET OpenMP::OpenMP_CXX AND OpenMP_FOUND)
    find_package (Threads REQUIRED)
    add_library (OpenMP::OpenMP_CXX IMPORTED INTERFACE)
    set_property (TARGET OpenMP::OpenMP_CXX PROPERTY INTERFACE_COMPILE_OPTIONS
                                                     ${OpenMP_CXX_FLAGS})
    set_property (
      TARGET OpenMP::OpenMP_CXX PROPERTY INTERFACE_LINK_LIBRARIES
                                         ${OpenMP_CXX_FLAGS} Threads::Threads)
  endif ()
endif ()

if (LIBCMAES_BUILD_PYTHON)
  find_package (PythonInterp)
  find_package (PythonLibs)
  if (NOT TARGET Python::Module AND PYTHONLIBS_FOUND)
    add_library (Python::Module SHARED IMPORTED)
    set_target_properties (
      Python::Module
      PROPERTIES IMPORTED_LOCATION ${PYTHON_LIBRARIES}
                 INTERFACE_INCLUDE_DIRECTORIES ${PYTHON_INCLUDE_DIRS})
  endif ()

  find_package (NumPy)

  # python site dir
  if (PYTHON_EXECUTABLE AND NOT DEFINED PYTHON_SITE_PACKAGES)
    # $ENV{VIRTUAL_ENV} does not work because not exported by activate :/
    execute_process (
      COMMAND
        ${PYTHON_EXECUTABLE} -c "
import sys;
if 'real_prefix' in dir(sys):
    print(sys.prefix)"
        # sys.real_prefix is only defined when in virtualenv
      OUTPUT_VARIABLE SITE_DIR_PREFIX
      RESULT_VARIABLE _exit_code
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    if (NOT ${SITE_DIR_PREFIX}) # if not set, back to previous behaviour
      set (SITE_DIR_PREFIX ${CMAKE_INSTALL_PREFIX})
    endif ()
    execute_process (
      COMMAND
        ${PYTHON_EXECUTABLE} -c
        "from distutils import sysconfig; print(sysconfig.get_python_lib(plat_specific=True, prefix='${SITE_DIR_PREFIX}'))"
      OUTPUT_VARIABLE _ABS_PYTHON_SITE_PACKAGES
      RESULT_VARIABLE _exit_code
      OUTPUT_STRIP_TRAILING_WHITESPACE)

    if (NOT ${_exit_code})
      get_filename_component (_ABS_PYTHON_SITE_PACKAGES
                              ${_ABS_PYTHON_SITE_PACKAGES} ABSOLUTE)
      file (RELATIVE_PATH _REL_PYTHON_SITE_PACKAGES ${CMAKE_INSTALL_PREFIX}
            ${_ABS_PYTHON_SITE_PACKAGES})
      set (PYTHON_SITE_PACKAGES ${_REL_PYTHON_SITE_PACKAGES})
    else ()
      message (SEND_ERROR "Could not run ${PYTHON_EXECUTABLE}")
    endif ()
  endif ()

  # boost-python
  find_package(Boost COMPONENTS python)
  if (Boost_FOUND)
    set (HAVE_BOOST_PYTHON TRUE)
    if (NOT TARGET Boost::python)
      add_library (Boost::python SHARED IMPORTED)
      set_target_properties (
        Boost::python
        PROPERTIES IMPORTED_LOCATION ${Boost_LIBRARIES}
                   INTERFACE_INCLUDE_DIRECTORIES ${Boost_INCLUDE_DIRS})
    endif ()
    message (STATUS "Found boost-python: ${Boost_LIBRARIES}")
  endif ()
endif ()

set (HAVE_SURROG ${LIBCMAES_ENABLE_SURROG})
configure_file (include/libcmaes/libcmaes_config.h.cmake.in include/libcmaes/libcmaes_config.h)
install (FILES ${PROJECT_BINARY_DIR}/include/libcmaes/libcmaes_config.h
  DESTINATION ${RELATIVE_INSTALL_INCLUDE_DIR}/libcmaes)

if (NOT MSVC)
  set (prefix ${CMAKE_INSTALL_PREFIX})
  set (exec_prefix ${prefix})
  set (libdir ${INSTALL_LIB_DIR})
  set (includedir ${INSTALL_INCLUDE_DIR})
  set (VERSION ${LIBCMAES_VERSION_STRING})
  configure_file (libcmaes.pc.in libcmaes.pc @ONLY)
  install (FILES ${CMAKE_CURRENT_BINARY_DIR}/libcmaes.pc
           DESTINATION ${RELATIVE_INSTALL_LIB_DIR}/pkgconfig)
endif ()

# ---------- subdirectories ----------
add_subdirectory (src)
if (LIBCMAES_BUILD_EXAMPLES)
  add_subdirectory (examples)
endif ()
if (LIBCMAES_BUILD_TESTS)
  enable_testing ()
  add_subdirectory (tests)
endif ()

if (Boost_FOUND AND NumPy_FOUND)
  add_subdirectory (python)
endif ()

# ------------ export ------------
add_library (libcmaes::cmaes ALIAS cmaes)

include (CMakePackageConfigHelpers)
write_basic_package_version_file (
  libcmaesConfigVersion.cmake VERSION ${PACKAGE_VERSION}
  COMPATIBILITY AnyNewerVersion)

# export target to build directory
export (TARGETS cmaes NAMESPACE libcmaes:: FILE libcmaesTargets.cmake)
# export target on install
install (
  EXPORT libcmaesTargets
  FILE libcmaesTargets.cmake
  NAMESPACE libcmaes::
  DESTINATION ${RELATIVE_INSTALL_CMAKE_DIR})

configure_file (libcmaesConfig.cmake.in libcmaesConfig.cmake @ONLY)
install (
  FILES "${CMAKE_CURRENT_BINARY_DIR}/libcmaesConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/libcmaesConfigVersion.cmake"
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindEigen3.cmake"
  DESTINATION ${RELATIVE_INSTALL_CMAKE_DIR})
