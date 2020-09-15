python_add_module (lcmaes lcmaes.cc)

target_link_libraries (lcmaes cmaes Python::Module Python::NumPy Boost::python)

if (APPLE)
  set_target_properties (lcmaes PROPERTIES LINK_FLAGS
                                           "-undefined dynamic_lookup")
endif ()

install (TARGETS lcmaes DESTINATION ${PYTHON_SITE_PACKAGES})

set (PYINSTALLCHECK_ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}")

if (LIBCMAES_BUILD_TESTS)
  macro (cmaes_add_pytest name)
    add_test (NAME ${name} COMMAND ${PYTHON_EXECUTABLE}
                                   ${CMAKE_CURRENT_SOURCE_DIR}/${name}.py)
    set_tests_properties (${name} PROPERTIES ENVIRONMENT
                                             "${PYINSTALLCHECK_ENVIRONMENT}")
  endmacro ()

  cmaes_add_pytest (ptest)
  cmaes_add_pytest (ptest_bounds)
  cmaes_add_pytest (ptests_ls)
endif ()
