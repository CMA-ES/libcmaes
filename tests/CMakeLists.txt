
include_directories (${PROJECT_SOURCE_DIR}/src)
include_directories (${PROJECT_BINARY_DIR}/src)

macro (cmaes_add_test name)
  add_executable (t_${name} ${name}.cc)
  target_link_libraries (t_${name} cmaes)
  add_test (NAME ${name} COMMAND t_${name})
  if (WIN32)
    set_tests_properties (${name} PROPERTIES ENVIRONMENT "PATH=${PROJECT_BINARY_DIR}\\src\\${CMAKE_BUILD_TYPE};$ENV{PATH}")  # to load dll
  endif ()
endmacro ()

cmaes_add_test (simple-test)
cmaes_add_test (edm)
