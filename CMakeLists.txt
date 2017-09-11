cmake_minimum_required (VERSION 2.8)

if (NOT DEFINED CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type")
endif ()

project (libcmaes)

list (APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)


option (BUILD_SHARED_LIBS "Build libcmaes as a shared library" ON)
option (BUILD_PYTHON "build python bindings" ON)
option (BUILD_TESTS "build tests" ON)
option (BUILD_EXAMPLES "build examples" ON)
option (USE_OPENMP "Use OpenMP for multithreading" ON)
option (ENABLE_SURROG "support for surrogates" ON)
set (HAVE_SURROG ${ENABLE_SURROG})
option (USE_COMPILE_FEATURES "use cmake>=3.1 cxx11 detection" ON)
option (LINK_PYTHON "link python libraries" ON)

# Offer the user the choice of overriding the installation directories
set (INSTALL_LIB_DIR     lib${LIB_SUFFIX} CACHE PATH "Installation directory for libraries")
set (INSTALL_BIN_DIR     bin              CACHE PATH "Installation directory for executables")
set (INSTALL_INCLUDE_DIR include          CACHE PATH "Installation directory for header files")
set (INSTALL_DATA_DIR    share/libcmaes       CACHE PATH "Installation directory for data files")
set (INSTALL_CMAKE_DIR   ${INSTALL_LIB_DIR}/cmake/libcmaes CACHE PATH "Installation directory for cmake config files")

# Make relative paths absolute (needed later on)
foreach (p LIB BIN INCLUDE DATA CMAKE)
  set (var INSTALL_${p}_DIR)
  set (RELATIVE_INSTALL_${p}_DIR ${INSTALL_${p}_DIR})
  if (NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif ()
endforeach ()


set (CMAKE_INSTALL_RPATH ${INSTALL_LIB_DIR})
set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
if (POLICY CMP0042)
  # Set MACOSX_RPATH to ON
  cmake_policy (SET CMP0042 NEW)
endif()

include (CheckIncludeFile)
include (CheckCXXCompilerFlag)

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

# eigen
find_package (Eigen3 NO_MODULE)
if (EIGEN3_FOUND)
  include (${EIGEN3_USE_FILE})
  message(STATUS "Eigen3 found: ${EIGEN3_INCLUDE_DIRS}")
else ()
  find_package (Eigen3 REQUIRED)
  include_directories (${EIGEN3_INCLUDE_DIRS})
endif()

# omp
if (USE_OPENMP)
  find_package (OpenMP)
  if (OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")
  endif ()
endif ()


if (CMAKE_VERSION VERSION_LESS 3.1 OR NOT USE_COMPILE_FEATURES)
  check_cxx_compiler_flag ("-std=c++11" SUPPORTS_STDCXX11)
  if (SUPPORTS_STDCXX11)
    set (CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")
  endif ()
endif ()

if (BUILD_PYTHON)
  find_package (PythonInterp)
  find_package (PythonLibs)
  find_package (NumPy)

  # python site dir
  if (PYTHON_EXECUTABLE AND NOT DEFINED PYTHON_SITE_PACKAGES)
    # $ENV{VIRTUAL_ENV} does not work because not exported by activate :/
    execute_process (COMMAND ${PYTHON_EXECUTABLE} -c "
import sys;
if 'real_prefix' in dir(sys):
    print(sys.prefix)"
    # sys.real_prefix is only defined when in virtualenv
      OUTPUT_VARIABLE SITE_DIR_PREFIX
      RESULT_VARIABLE _exit_code
      OUTPUT_STRIP_TRAILING_WHITESPACE )
    if (${SITE_DIR_PREFIX}) # if not set, back to previous behaviour
      set (SITE_DIR_PREFIX ${CMAKE_INSTALL_PREFIX})
    endif ()
    execute_process (COMMAND ${PYTHON_EXECUTABLE} -c "from distutils import sysconfig; print(sysconfig.get_python_lib(plat_specific=True, prefix='${SITE_DIR_PREFIX}'))"
                         OUTPUT_VARIABLE _ABS_PYTHON_SITE_PACKAGES
                         RESULT_VARIABLE _exit_code
                         OUTPUT_STRIP_TRAILING_WHITESPACE )

    if (NOT ${_exit_code})
      get_filename_component (_ABS_PYTHON_SITE_PACKAGES ${_ABS_PYTHON_SITE_PACKAGES} ABSOLUTE)
      file (RELATIVE_PATH _REL_PYTHON_SITE_PACKAGES ${CMAKE_INSTALL_PREFIX} ${_ABS_PYTHON_SITE_PACKAGES})
      set (PYTHON_SITE_PACKAGES ${_REL_PYTHON_SITE_PACKAGES})
    else ()
      message (SEND_ERROR "Could not run ${PYTHON_EXECUTABLE}")
    endif ()
  endif ()

  # boost-python
  set (BOOST_PYTHON_COMPONENT_NAMES
    python-py${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR} python${PYTHON_VERSION_MAJOR} python)
  foreach (NAME ${BOOST_PYTHON_COMPONENT_NAMES})
    if (NOT Boost_FOUND)
      message (STATUS "Looking for component ${NAME}")
      find_package(Boost QUIET COMPONENTS ${NAME})
    endif ()
  endforeach ()
  if (Boost_FOUND)
    set (HAVE_BOOST_PYTHON TRUE)
    message(STATUS "Found boost-python: ${Boost_LIBRARIES}")
  endif ()
endif ()

set (LIBCMAES_MAJOR_VERSION "0")
set (LIBCMAES_MINOR_VERSION "9")
set (LIBCMAES_BUGFIX_VERSION "5")
set (LIBCMAES_VERSION_STRING ${LIBCMAES_MAJOR_VERSION}.${LIBCMAES_MINOR_VERSION}.${LIBCMAES_BUGFIX_VERSION})



configure_file (libcmaes_config.h.cmake.in libcmaes_config.h)
install (FILES ${PROJECT_BINARY_DIR}/libcmaes_config.h DESTINATION ${RELATIVE_INSTALL_INCLUDE_DIR}/libcmaes)
include_directories (${CMAKE_CURRENT_BINARY_DIR})

if (NOT MSVC)
  set (prefix ${CMAKE_INSTALL_PREFIX})
  set (exec_prefix ${prefix})
  set (libdir ${INSTALL_LIB_DIR})
  set (includedir ${INSTALL_INCLUDE_DIR})
  set (VERSION ${LIBCMAES_VERSION_STRING})
  configure_file (libcmaes.pc.in libcmaes.pc @ONLY)
  install (FILES ${CMAKE_CURRENT_BINARY_DIR}/libcmaes.pc DESTINATION ${RELATIVE_INSTALL_LIB_DIR}/pkgconfig)
endif ()


enable_testing ()

add_subdirectory (src)
if (BUILD_EXAMPLES)
  add_subdirectory (examples)
endif ()
if (BUILD_TESTS)
  add_subdirectory (tests)
endif ()

if (Boost_FOUND AND NUMPY_FOUND)
  add_subdirectory (python)
endif ()
