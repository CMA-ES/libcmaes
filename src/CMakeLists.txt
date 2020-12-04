set (
  LIBCMAES_SOURCES
  acovarianceupdate.cc
  cmaparameters.cc
  cmastopcriteria.cc
  covarianceupdate.cc
  esostrategy.cc
  pwq_bound_strategy.cc
  vdcmaupdate.cc
  bipopcmastrategy.cc
  cmasolutions.cc
  cmastrategy.cc
  errstats.cc
  ipopcmastrategy.cc)

set(header_path "${PROJECT_SOURCE_DIR}/include/libcmaes")
set (LIBCMAES_HEADERS
  ${header_path}/cmaes.h
  ${header_path}/opti_err.h
  ${header_path}/eo_matrix.h
  ${header_path}/cmastrategy.h
  ${header_path}/esoptimizer.h
  ${header_path}/esostrategy.h
  ${header_path}/cmasolutions.h
  ${header_path}/parameters.h
  ${header_path}/cmaparameters.h
  ${header_path}/cmastopcriteria.h
  ${header_path}/ipopcmastrategy.h
  ${header_path}/bipopcmastrategy.h
  ${header_path}/covarianceupdate.h
  ${header_path}/acovarianceupdate.h
  ${header_path}/vdcmaupdate.h
  ${header_path}/pwq_bound_strategy.h
  ${header_path}/eigenmvn.h
  ${header_path}/candidate.h
  ${header_path}/genopheno.h
  ${header_path}/noboundstrategy.h
  ${header_path}/scaling.h
  ${header_path}/llogging.h
  ${header_path}/errstats.h
  ${header_path}/pli.h
  ${header_path}/contour.h)

if (LIBCMAES_ENABLE_SURROG)
  list (APPEND LIBCMAES_SOURCES surrogatestrategy.cc)
  list (APPEND LIBCMAES_HEADERS ${header_path}/surrcmaes.h ${header_path}/surrogatestrategy.h)
endif ()

if (LIBCMAES_BUILD_SHARED_LIBS)
  add_library (cmaes SHARED ${LIBCMAES_SOURCES})
else ()
  add_library (cmaes STATIC ${LIBCMAES_SOURCES})
endif ()

target_include_directories (
  cmaes PUBLIC $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
               $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  )
target_link_libraries (cmaes PUBLIC Eigen3::Eigen)
if (LIBCMAES_USE_OPENMP)
  target_link_libraries (cmaes PUBLIC OpenMP::OpenMP_CXX)
endif ()

target_compile_features (cmaes PUBLIC cxx_nonstatic_member_init)
if (${CMAKE_VERSION} VERSION_GREATER 3.8)
  target_compile_features (cmaes PUBLIC cxx_std_11)
  set_target_properties (cmaes PROPERTIES CXX_EXTENSIONS OFF)
else ()
  set_target_properties (cmaes PROPERTIES CXX_STANDARD 11 CXX_STANDARD_REQUIRED
                                          YES CXX_EXTENSIONS OFF)
endif ()

if (MSVC)
  target_compile_options (cmaes PRIVATE "/bigobj")
endif ()

set_target_properties (cmaes PROPERTIES SOVERSION 0 VERSION 0.0.0
                                        POSITION_INDEPENDENT_CODE ON)

include (GenerateExportHeader)
generate_export_header (cmaes EXPORT_FILE_NAME ${CMAKE_BINARY_DIR}/include/libcmaes/cmaes_export.h)

install (
  TARGETS cmaes
  EXPORT libcmaesTargets
  RUNTIME DESTINATION ${RELATIVE_INSTALL_BIN_DIR}
  LIBRARY DESTINATION ${RELATIVE_INSTALL_LIB_DIR}
  ARCHIVE DESTINATION ${RELATIVE_INSTALL_LIB_DIR}
  INCLUDES DESTINATION ${RELATIVE_INSTALL_INCLUDE_DIR})


install (FILES ${LIBCMAES_HEADERS}  ${CMAKE_BINARY_DIR}/include/libcmaes/cmaes_export.h
         DESTINATION ${RELATIVE_INSTALL_INCLUDE_DIR}/libcmaes)

if (LIBCMAES_ENABLE_SURROG)
  install (DIRECTORY ${header_path}/surrogates
           DESTINATION ${RELATIVE_INSTALL_INCLUDE_DIR}/libcmaes)
endif ()
