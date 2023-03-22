Python3_add_library (lcmaes MODULE lcmaes.cc)

if(MSVC)
  add_custom_command(TARGET lcmaes POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:cmaes> $<TARGET_FILE_DIR:lcmaes> )
endif()

target_link_libraries (lcmaes PRIVATE cmaes Python3::Module Python3::NumPy Boost::python)

if (APPLE)
  set_target_properties (lcmaes PROPERTIES LINK_FLAGS
                                           "-undefined dynamic_lookup")
endif ()

install (TARGETS lcmaes DESTINATION ${PYTHON_SITE_PACKAGES})
if(MSVC)
install(TARGETS cmaes DESTINATION ${PYTHON_SITE_PACKAGES})
endif()


if (LIBCMAES_BUILD_TESTS)
  macro (cmaes_add_pytest name)
    add_test (NAME ${name} COMMAND ${Python3_EXECUTABLE}
                                   ${CMAKE_CURRENT_SOURCE_DIR}/${name}.py)
    set_tests_properties (${name} PROPERTIES ENVIRONMENT
                                             "PYTHONPATH=$<TARGET_FILE_DIR:lcmaes>")
  endmacro ()

  cmaes_add_pytest (ptest)
  cmaes_add_pytest (ptest_bounds)
  cmaes_add_pytest (ptests_ls)
endif ()
