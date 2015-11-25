# Install script for directory: /Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Geometry

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
   "/Eigen/src/Geometry/AlignedBox.h;/Eigen/src/Geometry/AngleAxis.h;/Eigen/src/Geometry/EulerAngles.h;/Eigen/src/Geometry/Homogeneous.h;/Eigen/src/Geometry/Hyperplane.h;/Eigen/src/Geometry/OrthoMethods.h;/Eigen/src/Geometry/ParametrizedLine.h;/Eigen/src/Geometry/Quaternion.h;/Eigen/src/Geometry/Rotation2D.h;/Eigen/src/Geometry/RotationBase.h;/Eigen/src/Geometry/Scaling.h;/Eigen/src/Geometry/Transform.h;/Eigen/src/Geometry/Translation.h;/Eigen/src/Geometry/Umeyama.h")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/Eigen/src/Geometry" TYPE FILE FILES
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Geometry/AlignedBox.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Geometry/AngleAxis.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Geometry/EulerAngles.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Geometry/Homogeneous.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Geometry/Hyperplane.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Geometry/OrthoMethods.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Geometry/ParametrizedLine.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Geometry/Quaternion.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Geometry/Rotation2D.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Geometry/RotationBase.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Geometry/Scaling.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Geometry/Transform.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Geometry/Translation.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Geometry/Umeyama.h"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Geometry/arch/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/Geometry/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
