# Install script for directory: /Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/src/NonLinearOptimization

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
   "/unsupported/Eigen/src/NonLinearOptimization/chkder.h;/unsupported/Eigen/src/NonLinearOptimization/covar.h;/unsupported/Eigen/src/NonLinearOptimization/dogleg.h;/unsupported/Eigen/src/NonLinearOptimization/fdjac1.h;/unsupported/Eigen/src/NonLinearOptimization/HybridNonLinearSolver.h;/unsupported/Eigen/src/NonLinearOptimization/LevenbergMarquardt.h;/unsupported/Eigen/src/NonLinearOptimization/lmpar.h;/unsupported/Eigen/src/NonLinearOptimization/qrsolv.h;/unsupported/Eigen/src/NonLinearOptimization/r1mpyq.h;/unsupported/Eigen/src/NonLinearOptimization/r1updt.h;/unsupported/Eigen/src/NonLinearOptimization/rwupdt.h")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/unsupported/Eigen/src/NonLinearOptimization" TYPE FILE FILES
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/src/NonLinearOptimization/chkder.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/src/NonLinearOptimization/covar.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/src/NonLinearOptimization/dogleg.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/src/NonLinearOptimization/fdjac1.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/src/NonLinearOptimization/HybridNonLinearSolver.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/src/NonLinearOptimization/LevenbergMarquardt.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/src/NonLinearOptimization/lmpar.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/src/NonLinearOptimization/qrsolv.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/src/NonLinearOptimization/r1mpyq.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/src/NonLinearOptimization/r1updt.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/src/NonLinearOptimization/rwupdt.h"
    )
endif()

