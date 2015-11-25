# Install script for directory: /Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/src/SparseExtra

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
   "/unsupported/Eigen/src/SparseExtra/BlockOfDynamicSparseMatrix.h;/unsupported/Eigen/src/SparseExtra/BlockSparseMatrix.h;/unsupported/Eigen/src/SparseExtra/DynamicSparseMatrix.h;/unsupported/Eigen/src/SparseExtra/MarketIO.h;/unsupported/Eigen/src/SparseExtra/MatrixMarketIterator.h;/unsupported/Eigen/src/SparseExtra/RandomSetter.h")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/unsupported/Eigen/src/SparseExtra" TYPE FILE FILES
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/src/SparseExtra/BlockOfDynamicSparseMatrix.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/src/SparseExtra/BlockSparseMatrix.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/src/SparseExtra/DynamicSparseMatrix.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/src/SparseExtra/MarketIO.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/src/SparseExtra/MatrixMarketIterator.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/src/SparseExtra/RandomSetter.h"
    )
endif()

