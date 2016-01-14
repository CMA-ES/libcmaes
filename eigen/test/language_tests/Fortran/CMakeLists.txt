project(test NONE)
    cmake_minimum_required(VERSION 2.8.0)
    set (CMAKE_Fortran_FLAGS "")
    set (CMAKE_EXE_LINKER_FLAGS " -L/usr/local/opt/lapack/lib")
    enable_language(Fortran OPTIONAL)
  