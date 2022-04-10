# SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas Müller
#
# SPDX-License-Identifier: GPL-3.0-or-later

set(CMAKE_CXX_FLAGS_DEBUG
	  "${CMAKE_CXX_FLAGS_DEBUG} -Wall"
)
set(CMAKE_CXX_FLAGS_RELEASE
	  "${CMAKE_CXX_FLAGS_RELEASE} -Wall ")

if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS_DEBUG
       "${CMAKE_CXX_FLAGS_DEBUG} -Wextra")
  set(CMAKE_CXX_FLAGS_RELEASE
    "${CMAKE_CXX_FLAGS_RELEASE} -O3 -Wextra")
endif(CMAKE_COMPILER_IS_GNUCXX)

if (MSVC)
     set(CMAKE_CXX_FLAGS_DEBUG
       "${CMAKE_CXX_FLAGS_DEBUG} /permissive- /bigobj /w44101")

     set(CMAKE_CXX_FLAGS_RELEASE
        "${CMAKE_CXX_FLAGS_RELEASE} /permissive- /bigobj /w44101 /Ox")

     set(MSVC_DISABLED_WARNINGS_LIST
      "C4061" # enumerator 'identifier' in switch of enum 'enumeration' is not
              # explicitly handled by a case label
              # Disable this because it flags even when there is a default.
      "C4068"
      "C4100" # 'exarg' : unreferenced formal parameter
      "C4127" # conditional expression is constant
      "C4200" # nonstandard extension used : zero-sized array in
              # struct/union.
      "C4204" # nonstandard extension used: non-constant aggregate initializer
      "C4221" # nonstandard extension used : 'identifier' : cannot be
              # initialized using address of automatic variable
      "C4242" # 'function' : conversion from 'int' to 'uint8_t',
              # possible loss of data
      "C4244" # 'function' : conversion from 'int' to 'uint8_t',
              # possible loss of data
      "C4245" # 'initializing' : conversion from 'long' to
              # 'unsigned long', signed/unsigned mismatch
      "C4267" # conversion from 'size_t' to 'int', possible loss of data
      "C4355"
      "C4371" # layout of class may have changed from a previous version of the
              # compiler due to better packing of member '...'
      "C4388" # signed/unsigned mismatch
      "C4296" # '>=' : expression is always true
      "C4350" # behavior change: 'std::_Wrap_alloc...'
      "C4365" # '=' : conversion from 'size_t' to 'int',
              # signed/unsigned mismatch
      "C4389" # '!=' : signed/unsigned mismatch
      "C4464" # relative include path contains '..'
      "C4510" # 'argument' : default constructor could not be generated
      "C4571"
      "C4512" # 'argument' : assignment operator could not be generated
      "C4514" # 'function': unreferenced inline function has been removed
      "C4548" # expression before comma has no effect; expected expression with
              # side-effect" caused by FD_* macros.
      "C4610" # struct 'argument' can never be instantiated - user defined
              # constructor required.
      "C4619"
      "C4623" # default constructor was implicitly defined as deleted
      "C4625" # copy constructor could not be generated because a base class
              # copy constructor is inaccessible or deleted
      "C4626" # assignment operator could not be generated because a base class
              # assignment operator is inaccessible or deleted
    "C4643"
      "C4668" # 'symbol' is not defined as a preprocessor macro, replacing with
              # '0' for 'directives'
              # Disable this because GTest uses it everywhere.
      "C4706" # assignment within conditional expression
      "C4710" # 'function': function not inlined
      "C4711" # function 'function' selected for inline expansion
      "C4800" # 'int' : forcing value to bool 'true' or 'false'
              # (performance warning)
      "C4820" # 'bytes' bytes padding added after construct 'member_name'
    "C4868"
    "C4996"
      "C5026" # move constructor was implicitly defined as deleted
      "C5027" # move assignment operator was implicitly defined as deleted
      "C5031"
    "C5039"
      "C5045"
      )
      string(REPLACE "C" " -wd" MSVC_DISABLED_WARNINGS_STR
                            ${MSVC_DISABLED_WARNINGS_LIST})

      set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} ${MSVC_DISABLED_WARNINGS_STR}")
 endif()
