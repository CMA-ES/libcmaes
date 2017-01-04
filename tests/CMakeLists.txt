
include_directories (${PROJECT_SOURCE_DIR}/src)

macro (cmaes_add_test name)
  add_executable (t_${name} ${name}.cc)
  target_link_libraries (t_${name} cmaes)
  add_test (NAME ${name} COMMAND t_${name})
endmacro ()

cmaes_add_test (simple-test)
cmaes_add_test (edm)
