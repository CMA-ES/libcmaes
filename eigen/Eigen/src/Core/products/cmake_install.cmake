# Install script for directory: /Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/products

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
   "/Eigen/src/Core/products/GeneralBlockPanelKernel.h;/Eigen/src/Core/products/GeneralMatrixMatrix.h;/Eigen/src/Core/products/GeneralMatrixMatrix_MKL.h;/Eigen/src/Core/products/GeneralMatrixMatrixTriangular.h;/Eigen/src/Core/products/GeneralMatrixMatrixTriangular_MKL.h;/Eigen/src/Core/products/GeneralMatrixVector.h;/Eigen/src/Core/products/GeneralMatrixVector_MKL.h;/Eigen/src/Core/products/Parallelizer.h;/Eigen/src/Core/products/SelfadjointMatrixMatrix.h;/Eigen/src/Core/products/SelfadjointMatrixMatrix_MKL.h;/Eigen/src/Core/products/SelfadjointMatrixVector.h;/Eigen/src/Core/products/SelfadjointMatrixVector_MKL.h;/Eigen/src/Core/products/SelfadjointProduct.h;/Eigen/src/Core/products/SelfadjointRank2Update.h;/Eigen/src/Core/products/TriangularMatrixMatrix.h;/Eigen/src/Core/products/TriangularMatrixMatrix_MKL.h;/Eigen/src/Core/products/TriangularMatrixVector.h;/Eigen/src/Core/products/TriangularMatrixVector_MKL.h;/Eigen/src/Core/products/TriangularSolverMatrix.h;/Eigen/src/Core/products/TriangularSolverMatrix_MKL.h;/Eigen/src/Core/products/TriangularSolverVector.h")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/Eigen/src/Core/products" TYPE FILE FILES
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/products/GeneralBlockPanelKernel.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/products/GeneralMatrixMatrix.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/products/GeneralMatrixMatrix_MKL.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/products/GeneralMatrixMatrixTriangular.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/products/GeneralMatrixMatrixTriangular_MKL.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/products/GeneralMatrixVector.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/products/GeneralMatrixVector_MKL.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/products/Parallelizer.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/products/SelfadjointMatrixMatrix.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/products/SelfadjointMatrixMatrix_MKL.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/products/SelfadjointMatrixVector.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/products/SelfadjointMatrixVector_MKL.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/products/SelfadjointProduct.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/products/SelfadjointRank2Update.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/products/TriangularMatrixMatrix.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/products/TriangularMatrixMatrix_MKL.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/products/TriangularMatrixVector.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/products/TriangularMatrixVector_MKL.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/products/TriangularSolverMatrix.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/products/TriangularSolverMatrix_MKL.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/products/TriangularSolverVector.h"
    )
endif()

