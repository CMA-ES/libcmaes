# Install script for directory: /Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/src/IterativeSolvers

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
   "/unsupported/Eigen/src/IterativeSolvers/ConstrainedConjGrad.h;/unsupported/Eigen/src/IterativeSolvers/DGMRES.h;/unsupported/Eigen/src/IterativeSolvers/GMRES.h;/unsupported/Eigen/src/IterativeSolvers/IncompleteLU.h;/unsupported/Eigen/src/IterativeSolvers/IterationController.h;/unsupported/Eigen/src/IterativeSolvers/MINRES.h;/unsupported/Eigen/src/IterativeSolvers/Scaling.h")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/unsupported/Eigen/src/IterativeSolvers" TYPE FILE FILES
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/src/IterativeSolvers/ConstrainedConjGrad.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/src/IterativeSolvers/DGMRES.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/src/IterativeSolvers/GMRES.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/src/IterativeSolvers/IncompleteLU.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/src/IterativeSolvers/IterationController.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/src/IterativeSolvers/MINRES.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/src/IterativeSolvers/Scaling.h"
    )
endif()

