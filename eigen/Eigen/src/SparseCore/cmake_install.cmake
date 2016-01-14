# Install script for directory: /Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseCore

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
   "/Eigen/src/SparseCore/AmbiVector.h;/Eigen/src/SparseCore/CompressedStorage.h;/Eigen/src/SparseCore/ConservativeSparseSparseProduct.h;/Eigen/src/SparseCore/MappedSparseMatrix.h;/Eigen/src/SparseCore/SparseAssign.h;/Eigen/src/SparseCore/SparseBlock.h;/Eigen/src/SparseCore/SparseColEtree.h;/Eigen/src/SparseCore/SparseCompressedBase.h;/Eigen/src/SparseCore/SparseCwiseBinaryOp.h;/Eigen/src/SparseCore/SparseCwiseUnaryOp.h;/Eigen/src/SparseCore/SparseDenseProduct.h;/Eigen/src/SparseCore/SparseDiagonalProduct.h;/Eigen/src/SparseCore/SparseDot.h;/Eigen/src/SparseCore/SparseFuzzy.h;/Eigen/src/SparseCore/SparseMap.h;/Eigen/src/SparseCore/SparseMatrix.h;/Eigen/src/SparseCore/SparseMatrixBase.h;/Eigen/src/SparseCore/SparsePermutation.h;/Eigen/src/SparseCore/SparseProduct.h;/Eigen/src/SparseCore/SparseRedux.h;/Eigen/src/SparseCore/SparseRef.h;/Eigen/src/SparseCore/SparseSelfAdjointView.h;/Eigen/src/SparseCore/SparseSolverBase.h;/Eigen/src/SparseCore/SparseSparseProductWithPruning.h;/Eigen/src/SparseCore/SparseTranspose.h;/Eigen/src/SparseCore/SparseTriangularView.h;/Eigen/src/SparseCore/SparseUtil.h;/Eigen/src/SparseCore/SparseVector.h;/Eigen/src/SparseCore/SparseView.h;/Eigen/src/SparseCore/TriangularSolver.h")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/Eigen/src/SparseCore" TYPE FILE FILES
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseCore/AmbiVector.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseCore/CompressedStorage.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseCore/ConservativeSparseSparseProduct.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseCore/MappedSparseMatrix.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseCore/SparseAssign.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseCore/SparseBlock.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseCore/SparseColEtree.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseCore/SparseCompressedBase.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseCore/SparseCwiseBinaryOp.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseCore/SparseCwiseUnaryOp.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseCore/SparseDenseProduct.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseCore/SparseDiagonalProduct.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseCore/SparseDot.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseCore/SparseFuzzy.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseCore/SparseMap.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseCore/SparseMatrix.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseCore/SparseMatrixBase.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseCore/SparsePermutation.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseCore/SparseProduct.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseCore/SparseRedux.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseCore/SparseRef.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseCore/SparseSelfAdjointView.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseCore/SparseSolverBase.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseCore/SparseSparseProductWithPruning.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseCore/SparseTranspose.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseCore/SparseTriangularView.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseCore/SparseUtil.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseCore/SparseVector.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseCore/SparseView.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseCore/TriangularSolver.h"
    )
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseCore/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
