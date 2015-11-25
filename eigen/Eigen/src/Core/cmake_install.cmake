# Install script for directory: /Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core

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
   "/Eigen/src/Core/Array.h;/Eigen/src/Core/ArrayBase.h;/Eigen/src/Core/ArrayWrapper.h;/Eigen/src/Core/Assign.h;/Eigen/src/Core/Assign_MKL.h;/Eigen/src/Core/AssignEvaluator.h;/Eigen/src/Core/BandMatrix.h;/Eigen/src/Core/Block.h;/Eigen/src/Core/BooleanRedux.h;/Eigen/src/Core/CommaInitializer.h;/Eigen/src/Core/CoreEvaluators.h;/Eigen/src/Core/CoreIterators.h;/Eigen/src/Core/CwiseBinaryOp.h;/Eigen/src/Core/CwiseNullaryOp.h;/Eigen/src/Core/CwiseUnaryOp.h;/Eigen/src/Core/CwiseUnaryView.h;/Eigen/src/Core/DenseBase.h;/Eigen/src/Core/DenseCoeffsBase.h;/Eigen/src/Core/DenseStorage.h;/Eigen/src/Core/Diagonal.h;/Eigen/src/Core/DiagonalMatrix.h;/Eigen/src/Core/DiagonalProduct.h;/Eigen/src/Core/Dot.h;/Eigen/src/Core/EigenBase.h;/Eigen/src/Core/ForceAlignedAccess.h;/Eigen/src/Core/Fuzzy.h;/Eigen/src/Core/GeneralProduct.h;/Eigen/src/Core/GenericPacketMath.h;/Eigen/src/Core/GlobalFunctions.h;/Eigen/src/Core/Inverse.h;/Eigen/src/Core/IO.h;/Eigen/src/Core/Map.h;/Eigen/src/Core/MapBase.h;/Eigen/src/Core/MathFunctions.h;/Eigen/src/Core/Matrix.h;/Eigen/src/Core/MatrixBase.h;/Eigen/src/Core/NestByValue.h;/Eigen/src/Core/NoAlias.h;/Eigen/src/Core/NumTraits.h;/Eigen/src/Core/PermutationMatrix.h;/Eigen/src/Core/PlainObjectBase.h;/Eigen/src/Core/Product.h;/Eigen/src/Core/ProductEvaluators.h;/Eigen/src/Core/Random.h;/Eigen/src/Core/Redux.h;/Eigen/src/Core/Ref.h;/Eigen/src/Core/Replicate.h;/Eigen/src/Core/ReturnByValue.h;/Eigen/src/Core/Reverse.h;/Eigen/src/Core/Select.h;/Eigen/src/Core/SelfAdjointView.h;/Eigen/src/Core/SelfCwiseBinaryOp.h;/Eigen/src/Core/Solve.h;/Eigen/src/Core/SolveTriangular.h;/Eigen/src/Core/StableNorm.h;/Eigen/src/Core/Stride.h;/Eigen/src/Core/Swap.h;/Eigen/src/Core/Transpose.h;/Eigen/src/Core/Transpositions.h;/Eigen/src/Core/TriangularMatrix.h;/Eigen/src/Core/VectorBlock.h;/Eigen/src/Core/VectorwiseOp.h;/Eigen/src/Core/Visitor.h")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/Eigen/src/Core" TYPE FILE FILES
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/Array.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/ArrayBase.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/ArrayWrapper.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/Assign.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/Assign_MKL.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/AssignEvaluator.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/BandMatrix.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/Block.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/BooleanRedux.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/CommaInitializer.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/CoreEvaluators.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/CoreIterators.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/CwiseBinaryOp.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/CwiseNullaryOp.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/CwiseUnaryOp.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/CwiseUnaryView.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/DenseBase.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/DenseCoeffsBase.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/DenseStorage.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/Diagonal.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/DiagonalMatrix.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/DiagonalProduct.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/Dot.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/EigenBase.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/ForceAlignedAccess.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/Fuzzy.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/GeneralProduct.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/GenericPacketMath.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/GlobalFunctions.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/Inverse.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/IO.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/Map.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/MapBase.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/MathFunctions.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/Matrix.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/MatrixBase.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/NestByValue.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/NoAlias.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/NumTraits.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/PermutationMatrix.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/PlainObjectBase.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/Product.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/ProductEvaluators.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/Random.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/Redux.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/Ref.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/Replicate.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/ReturnByValue.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/Reverse.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/Select.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/SelfAdjointView.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/SelfCwiseBinaryOp.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/Solve.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/SolveTriangular.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/StableNorm.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/Stride.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/Swap.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/Transpose.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/Transpositions.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/TriangularMatrix.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/VectorBlock.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/VectorwiseOp.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/Visitor.h"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/products/cmake_install.cmake")
  include("/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/util/cmake_install.cmake")
  include("/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/arch/cmake_install.cmake")
  include("/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/functors/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Core/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
