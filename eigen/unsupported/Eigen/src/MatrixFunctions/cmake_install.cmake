# Install script for directory: /Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/src/MatrixFunctions

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
   "/unsupported/Eigen/src/MatrixFunctions/MatrixExponential.h;/unsupported/Eigen/src/MatrixFunctions/MatrixFunction.h;/unsupported/Eigen/src/MatrixFunctions/MatrixLogarithm.h;/unsupported/Eigen/src/MatrixFunctions/MatrixPower.h;/unsupported/Eigen/src/MatrixFunctions/MatrixSquareRoot.h;/unsupported/Eigen/src/MatrixFunctions/StemFunction.h")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/unsupported/Eigen/src/MatrixFunctions" TYPE FILE FILES
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/src/MatrixFunctions/MatrixExponential.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/src/MatrixFunctions/MatrixFunction.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/src/MatrixFunctions/MatrixLogarithm.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/src/MatrixFunctions/MatrixPower.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/src/MatrixFunctions/MatrixSquareRoot.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/src/MatrixFunctions/StemFunction.h"
    )
endif()

