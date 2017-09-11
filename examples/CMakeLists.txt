
include_directories (${PROJECT_SOURCE_DIR}/src)
include_directories (${PROJECT_BINARY_DIR}/src)

macro (cmaes_add_example name)
  add_executable (sample-${name} sample-${name}.cc)
  target_link_libraries (sample-${name} cmaes)
endmacro ()

cmaes_add_example (code)
cmaes_add_example (code-genopheno)
cmaes_add_example (code-pfunc)
cmaes_add_example (code-ask-tell)
cmaes_add_example (code-ask-tell-uh)
cmaes_add_example (code-bounds)
cmaes_add_example (code-gradient)
cmaes_add_example (code-lscaling)
cmaes_add_example (code-lscaling-sigmas)
cmaes_add_example (code-pffunc)
