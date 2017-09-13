
include_directories (${PROJECT_SOURCE_DIR}/src)
include_directories (${PROJECT_BINARY_DIR}/src)
include_directories (${PYTHON_INCLUDE_DIR})
include_directories (${NUMPY_INCLUDE_DIRS})
include_directories (${Boost_INCLUDE_DIRS})

python_add_module (lcmaes lcmaes.cc)
target_link_libraries(lcmaes cmaes ${Boost_LIBRARIES})
if (LINK_PYTHON)
  target_link_libraries(lcmaes ${PYTHON_LIBRARIES})
endif ()

install (TARGETS lcmaes DESTINATION ${PYTHON_SITE_PACKAGES})

set (PYINSTALLCHECK_ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}")


macro (cmaes_add_pytest name)
  add_test (NAME ${name} COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/${name}.py)
  set_tests_properties (${name} PROPERTIES ENVIRONMENT "${PYINSTALLCHECK_ENVIRONMENT}")
endmacro ()

cmaes_add_pytest (ptest)
cmaes_add_pytest (ptest_bounds)
cmaes_add_pytest (ptests_ls)
