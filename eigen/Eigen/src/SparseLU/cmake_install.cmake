# Install script for directory: /Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseLU

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
   "/Eigen/src/SparseLU/SparseLU.h;/Eigen/src/SparseLU/SparseLU_column_bmod.h;/Eigen/src/SparseLU/SparseLU_column_dfs.h;/Eigen/src/SparseLU/SparseLU_copy_to_ucol.h;/Eigen/src/SparseLU/SparseLU_gemm_kernel.h;/Eigen/src/SparseLU/SparseLU_heap_relax_snode.h;/Eigen/src/SparseLU/SparseLU_kernel_bmod.h;/Eigen/src/SparseLU/SparseLU_Memory.h;/Eigen/src/SparseLU/SparseLU_panel_bmod.h;/Eigen/src/SparseLU/SparseLU_panel_dfs.h;/Eigen/src/SparseLU/SparseLU_pivotL.h;/Eigen/src/SparseLU/SparseLU_pruneL.h;/Eigen/src/SparseLU/SparseLU_relax_snode.h;/Eigen/src/SparseLU/SparseLU_Structs.h;/Eigen/src/SparseLU/SparseLU_SupernodalMatrix.h;/Eigen/src/SparseLU/SparseLU_Utils.h;/Eigen/src/SparseLU/SparseLUImpl.h")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/Eigen/src/SparseLU" TYPE FILE FILES
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseLU/SparseLU.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseLU/SparseLU_column_bmod.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseLU/SparseLU_column_dfs.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseLU/SparseLU_copy_to_ucol.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseLU/SparseLU_gemm_kernel.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseLU/SparseLU_heap_relax_snode.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseLU/SparseLU_kernel_bmod.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseLU/SparseLU_Memory.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseLU/SparseLU_panel_bmod.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseLU/SparseLU_panel_dfs.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseLU/SparseLU_pivotL.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseLU/SparseLU_pruneL.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseLU/SparseLU_relax_snode.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseLU/SparseLU_Structs.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseLU/SparseLU_SupernodalMatrix.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseLU/SparseLU_Utils.h"
    "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseLU/SparseLUImpl.h"
    )
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/Users/stevenjames/Documents/libcmaes/eigen/Eigen/src/SparseLU/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
