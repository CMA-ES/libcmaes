macro (cmaes_add_test name)
  add_executable (${name} ${name}.cc)
  target_link_libraries (${name} cmaes gflags)
  add_test (NAME ${name} COMMAND t_${name})
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
