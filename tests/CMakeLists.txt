find_package(gflags REQUIRED)

if (TARGET gflags::gflags)
  set (LIBCMAES_GFLAGS_TARGET gflags::gflags)
elseif (TARGET gflags)
  set (LIBCMAES_GFLAGS_TARGET gflags)
elseif (DEFINED GFLAGS_TARGET)
  set (LIBCMAES_GFLAGS_TARGET ${GFLAGS_TARGET})
else ()
  message (FATAL_ERROR "Could not determine a usable gflags target.")
endif ()

macro (cmaes_add_test name)
  add_executable (${name} ${name}.cc)
  target_link_libraries (${name} cmaes ${LIBCMAES_GFLAGS_TARGET})
  if(MSVC)
    target_compile_definitions(${name} PUBLIC _USE_MATH_DEFINES)
  endif()
  add_test (NAME ${name} COMMAND ${name})
  if (WIN32)
    set_tests_properties (
      ${name}
      PROPERTIES
        ENVIRONMENT
        "PATH=${PROJECT_BINARY_DIR}\\src\\${CMAKE_BUILD_TYPE}\;${PROJECT_BINARY_DIR}\\src\;$ENV{PATH}"
    )
  endif ()
endmacro ()

cmaes_add_test (simple-test)
cmaes_add_test (edm)
cmaes_add_test (test-functions)
