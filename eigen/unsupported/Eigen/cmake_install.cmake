# Install script for directory: /Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen

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
   "/unsupported/Eigen/AdolcForward;/unsupported/Eigen/AlignedVector3;/unsupported/Eigen/ArpackSupport;/unsupported/Eigen/AutoDiff;/unsupported/Eigen/BVH;/unsupported/Eigen/FFT;/unsupported/Eigen/IterativeSolvers;/unsupported/Eigen/KroneckerProduct;/unsupported/Eigen/LevenbergMarquardt;/unsupported/Eigen/MatrixFunctions;/unsupported/Eigen/MoreVectorization;/unsupported/Eigen/MPRealSupport;/unsupported/Eigen/NonLinearOptimization;/unsupported/Eigen/NumericalDiff;/unsupported/Eigen/OpenGLSupport;/unsupported/Eigen/Polynomials;/unsupported/Eigen/Skyline;/unsupported/Eigen/SparseExtra;/unsupported/Eigen/Splines")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/unsupported/Eigen" TYPE FILE FILES
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/AdolcForward"
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/AlignedVector3"
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/ArpackSupport"
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/AutoDiff"
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/BVH"
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/FFT"
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/IterativeSolvers"
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/KroneckerProduct"
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/LevenbergMarquardt"
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/MatrixFunctions"
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/MoreVectorization"
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/MPRealSupport"
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/NonLinearOptimization"
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/NumericalDiff"
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/OpenGLSupport"
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/Polynomials"
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/Skyline"
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/SparseExtra"
    "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/Splines"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/src/cmake_install.cmake")
  include("/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/CXX11/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/Users/stevenjames/Documents/libcmaes/eigen/unsupported/Eigen/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
