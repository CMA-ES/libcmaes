# - Try to find eigen lib
# Once done this will define
#
#  EIGEN3_FOUND - system has eigen lib
#  EIGEN3_INCLUDE_DIRS - the eigen include directories
#
#=============================================================================
# Copyright 2014, Julien Schueller
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distributed this file outside of CMake, substitute the full
#  License text for the above reference.)

find_path (EIGEN3_INCLUDE_DIR
           NAMES Eigen/Core
           PATH_SUFFIXES eigen3
          )

set (EIGEN3_INCLUDE_DIRS ${EIGEN3_INCLUDE_DIR})

# handle the QUIETLY and REQUIRED arguments
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (Eigen3 DEFAULT_MSG EIGEN3_INCLUDE_DIR)

mark_as_advanced (
  EIGEN3_INCLUDE_DIR
  EIGEN3_INCLUDE_DIRS
)
