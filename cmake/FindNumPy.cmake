# This module finds NumPy if it is installed, and sets the following variables
# indicating where it is.
#
# * NUMPY_FOUND         - was NumPy found
# * NUMPY_VERSION       - the version of NumPy found as a string
# * NUMPY_INCLUDE_DIRS  - path to the NumPy include files
#
# Additionally the imported target
#
# * Python::NumPy
#
# is defined.

# Finding NumPy involves calling the Python interpreter
if (NumPy_FIND_REQUIRED)
  find_package (PythonInterp REQUIRED)
else ()
  find_package (PythonInterp ${NumPy_FIND_REQUIRED})
endif ()

if (PYTHONINTERP_FOUND)
  execute_process (
    COMMAND ${PYTHON_EXECUTABLE} -c "import numpy as n; print(n.__version__)"
    RESULT_VARIABLE NumPy_VERSION_ERROR
    OUTPUT_VARIABLE NumPy_VERSION_STRING
    ERROR_VARIABLE NumPy_VERSION_ERROR_OUTPUT
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  execute_process (
    COMMAND ${PYTHON_EXECUTABLE} -c "import numpy as n; print(n.get_include())"
    RESULT_VARIABLE NumPy_INCLUDE_ERROR
    OUTPUT_VARIABLE NumPy_INCLUDE_DIRS
    ERROR_VARIABLE NumPy_INCLUDE_ERROR_OUTPUT
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  if (NOT NumPy_VERSION_ERROR)
    string (REGEX MATCH "^[0-9]+\\.[0-9]+\\.[0-9]+" NumPy_VERSION
                  ${NumPy_VERSION_STRING}
    )
    if ("${NumPy_VERSION}" STREQUAL "")
      set (NumPy_VERSION_ERROR YES)
    endif ()
  endif ()

  if (NOT NumPy_INCLUDE_ERROR) # canonicalize paths
    get_filename_component (NumPy_INCLUDE_DIRS ${NumPy_INCLUDE_DIRS} ABSOLUTE)
  endif ()

  if (NOT NumPy_VERSION_ERROR AND NOT NumPy_INCLUDE_ERROR)
    set (NumPy_FIND_SUCCESS TRUE)
  endif ()
endif ()

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (
  NumPy
  FOUND_VAR
  NumPy_FOUND
  REQUIRED_VARS
  NumPy_INCLUDE_DIRS
  NumPy_FIND_SUCCESS
  VERSION_VAR
  NumPy_VERSION
)

if (NumPy_FOUND AND NOT TARGET Python::NumPy)
  add_library (Python::NumPy INTERFACE IMPORTED)
  target_include_directories (Python::NumPy INTERFACE ${NumPy_INCLUDE_DIRS})
endif ()
