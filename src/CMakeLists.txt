
set (LIBCMAES_SOURCES
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
ipopcmastrategy.cc
surrogatestrategy.cc)

set (LIBCMAES_HEADERS cmaes.h opti_err.h eo_matrix.h cmastrategy.h esoptimizer.h esostrategy.h cmasolutions.h parameters.h cmaparameters.h cmastopcriteria.h ipopcmastrategy.h bipopcmastrategy.h covarianceupdate.h acovarianceupdate.h vdcmaupdate.h pwq_bound_strategy.h eigenmvn.h candidate.h genopheno.h noboundstrategy.h scaling.h llogging.h errstats.h pli.h contour.h)

if (ENABLE_SURROG)
  list (APPEND LIBCMAES_SOURCES surrogatestrategy.cc)
  list (APPEND LIBCMAES_HEADERS surrcmaes.h surrogatestrategy.h)
endif ()

if(LIBCMAES_BUILD_SHARED_LIBS)
  add_library (cmaes SHARED ${LIBCMAES_SOURCES})
else()
  add_library (cmaes STATIC ${LIBCMAES_SOURCES})
endif()

target_include_directories(cmaes PUBLIC ${PROJECT_BINARY_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(cmaes PUBLIC Eigen3::Eigen)

target_compile_features (cmaes PUBLIC cxx_nonstatic_member_init)
if (${CMAKE_VERSION} VERSION_GREATER 3.8)
  target_compile_features(cmaes PUBLIC cxx_std_11)
  set_target_properties(cmaes PROPERTIES CXX_EXTENSIONS OFF)
else ()
  set_target_properties(cmaes PROPERTIES
      CXX_STANDARD 11
      CXX_STANDARD_REQUIRED YES
      CXX_EXTENSIONS OFF
  )
endif ()

if(MSVC)
  target_compile_options(cmaes PRIVATE "/bigobj")
endif()

set_target_properties (cmaes PROPERTIES SOVERSION 0 VERSION 0.0.0 POSITION_INDEPENDENT_CODE ON)

include (GenerateExportHeader)
generate_export_header (cmaes)

install ( TARGETS cmaes
          EXPORT  LibcmaesLibraryDepends
          RUNTIME DESTINATION ${RELATIVE_INSTALL_BIN_DIR}
          LIBRARY DESTINATION ${RELATIVE_INSTALL_LIB_DIR}
          ARCHIVE DESTINATION ${RELATIVE_INSTALL_LIB_DIR}
        )

install (FILES ${LIBCMAES_HEADERS} ${CMAKE_CURRENT_BINARY_DIR}/cmaes_export.h
         DESTINATION ${RELATIVE_INSTALL_INCLUDE_DIR}/libcmaes)

if (ENABLE_SURROG)
  install (DIRECTORY surrogates DESTINATION ${RELATIVE_INSTALL_INCLUDE_DIR}/libcmaes)
endif ()
