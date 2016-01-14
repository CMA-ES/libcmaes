# Install script for directory: /Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/util

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Devel")
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/Eigen/src/Core/util/BlasUtil.h;/Eigen/src/Core/util/Constants.h;/Eigen/src/Core/util/DisableStupidWarnings.h;/Eigen/src/Core/util/ForwardDeclarations.h;/Eigen/src/Core/util/Macros.h;/Eigen/src/Core/util/Memory.h;/Eigen/src/Core/util/Meta.h;/Eigen/src/Core/util/MKL_support.h;/Eigen/src/Core/util/NonMPL2.h;/Eigen/src/Core/util/ReenableStupidWarnings.h;/Eigen/src/Core/util/StaticAssert.h;/Eigen/src/Core/util/XprHelper.h")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/Eigen/src/Core/util" TYPE FILE FILES
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/util/BlasUtil.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/util/Constants.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/util/DisableStupidWarnings.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/util/ForwardDeclarations.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/util/Macros.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/util/Memory.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/util/Meta.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/util/MKL_support.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/util/NonMPL2.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/util/ReenableStupidWarnings.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/util/StaticAssert.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/util/XprHelper.h"
    )
endif()

