# configure for libcmaes
# author
AC_INIT(libcmaes,0.9.5)
AC_MSG_NOTICE([libcmaes $PACKAGE_VERSION])
AC_CONFIG_SRCDIR(Makefile.am)
AM_INIT_AUTOMAKE([subdir-objects])
AM_CONFIG_HEADER([libcmaes_config.h])
AC_CONFIG_FILES(libcmaes.pc)

#ac_default_prefix=`pwd`
#bindir=$ac_default_prefix/bin/$host
#libdir=$ac_default_prefix/lib/$host

AC_PROG_INSTALL
AC_CANONICAL_HOST

AC_PROG_CC(clang++ g++ c++) # force c++ compiler on C files, otherwise BBOB lib doesn't link properly
AC_PROG_CXX(clang++ g++ c++)
AM_CONDITIONAL(HAVE_CLANG, test x"$CXX" = xclang++)
AC_PROG_INSTALL
AC_PROG_LIBTOOL

m4_include([m4/ax_cxx_compile_stdcxx_11.m4])
AX_CXX_COMPILE_STDCXX_11()

AC_LANG_PUSH([C++])
AC_MSG_CHECKING([for eigen3])
AC_ARG_WITH(eigen3-include,
AC_HELP_STRING([--with-eigen3-include=DIR], [Path to the eigen3 header files]),
					    [eig_incs=$withval],
					    [eig_incs=/usr/include/eigen3])
AC_CHECK_FILE($eig_incs/Eigen,,
	[AC_CHECK_FILE(/usr/local/include/eigen3,[eig_incs=/usr/local/include/eigen3])])

eig_ac_save_cppflags=$CPPFLAGS
CPPFLAGS="$CPPFLAGS -I$eig_incs"

lib_eig_avail=no
AC_CHECK_HEADER(Eigen/Dense,
	lib_eig_avail=yes,)

CPPFLAGS=$eig_ac_save_cppflags

if test x"$lib_eig_avail" = xyes; then
   AC_SUBST([EIGEN3_INC], [$eig_incs])
else
AC_MSG_ERROR([Cannot find Eigen3])
fi
AM_CONDITIONAL(HAVE_EIGEN3, test x"$lib_eig_avail" = xyes)

# Check if there is google-gflags library installed.
AC_ARG_WITH(gflags, AS_HELP_STRING[--with-gflags=GFLAGS_DIR],
  GFLAGS_CFLAGS="-I${with_gflags}/include"
  GFLAGS_LIBS="-L${with_gflags}/lib -lgflags"
  CFLAGS="$CFLAGS $GFLAGS_CFLAGS"
  LIBS="$LIBS $GFLAGS_LIBS"
)
AC_CHECK_LIB(gflags, main, ac_cv_have_libgflags=1, ac_cv_have_libgflags=0)
if test x"$ac_cv_have_libgflags" = x"1"; then
   AC_DEFINE(HAVE_LIB_GFLAGS, 1, [define if you have google gflags library])
   if test x"$GFLAGS_LIBS" = x""; then
      GFLAGS_LIBS="-lgflags"
   fi
else
  GFLAGS_CFLAGS=
  GFLAGS_LIBS=
fi
AC_SUBST(GFLAGS_LIBS)
AC_SUBST(GFLAGS_CFLAGS)
AM_CONDITIONAL(HAVE_GFLAGS, test x"$ac_cv_have_libgflags" = x"1")

# Check if there is google-glog library installed.
AC_ARG_ENABLE(gglog,
AC_HELP_STRING([--enable-gglog],[Enable google glog (default is NO)]),
[case "${enableval}" in 
	yes) gglog=true ;;
	no) gglog=false ;;
	*) AC_MSG_ERROR(bad_value ${enableval} for --enable-gglog) ;;
esac],[gglog=false])
if test x"$gglog" = x"true"; then
   AC_DEFINE(HAVE_GLOG,1,[google glog requested])
fi

AC_ARG_WITH(glog, AS_HELP_STRING[--with-glog=GLOG_DIR],
  GLOG_CFLAGS="-I${with_glog}/include"
  GLOG_LIBS="-L${with_glog}/lib -lglog"
  CFLAGS="$CFLAGS $GLOG_CFLAGS"
  LIBS="$LIBS $GLOG_LIBS"
)
AC_CHECK_LIB(glog, main, ac_cv_have_libglog=1, ac_cv_have_libglog=0)
if test x"$ac_cv_have_libglog" = x"1"; then
   AC_DEFINE(HAVE_LIB_GLOG, 1, [define if you have google glog library])
   if test x"$GLOG_LIBS" = x""; then
      GLOG_LIBS="-lglog"
   fi
else
  GLOG_CFLAGS=
  GLOG_LIBS=
fi
AC_SUBST(GLOG_LIBS)

# Check for gtest.
AC_CHECK_HEADER([gtest/gtest.h], gtest_header_avail=yes,)
AM_CONDITIONAL(HAVE_GTEST, test x"$gtest_header_avail" = xyes)

# Whether to compile with BBOB
AC_ARG_ENABLE(bbob,
AC_HELP_STRING([--enable-bbob],[Enable support for BBOB'2013 benchmark (default is NO)]),
[case "${enableval}" in 
	yes) bbob=true ;;
	no) bbob=false ;;
	*) AC_MSG_ERROR(bad_value ${enableval} for --enable-bbob) ;;
esac],[bbob=false])
AM_CONDITIONAL(HAVE_BBOB, test x$bbob = xtrue)

# Whether to compile with python wrapper support
AC_ARG_ENABLE(python,
AC_HELP_STRING([--enable-python],[Enable support for python wrapper (default is NO)]),
[case "${enableval}" in 
	yes) python=true ;;
	no) python=false ;;
	*) AC_MSG_ERROR(bad_value ${enableval} for --enable-python) ;;
esac],[python=false])
AM_CONDITIONAL(HAVE_PYTHON, test x$python = xtrue)
if test x"$python" = x"true"; then
#   AC_DEFINE(HAVE_PYTHON, 1, [python boost])
    AX_BOOST_PYTHON()
    if test "$BOOST_PYTHON_LIB" == ""; then
       AC_MSG_ERROR(boost python required for building the python wrapper,1)
    fi
fi

# Check for numpy
if test x"$python" = x"true"; then
AC_MSG_CHECKING([for Numpy include directory])
CPPFLAGS=$PYTHON_CPPFLAGS $CPPFLAGS
NUMPY_INCLUDE_DIR=`echo "import numpy; print numpy.get_include()" | $PYTHON - 2>/dev/null`
AC_SUBST(NUMPY_INCLUDE_DIR)
AC_CHECK_HEADER([${NUMPY_INCLUDE_DIR}/numpy/arrayobject.h],
                [NUMPY_HEADER=yes],
			[AC_MSG_ERROR([Numpy extension header not found])],
			   		          [#include "Python.h"])
AC_SUBST(NUMPY_HEADER)
fi

# Whether to only compile the library, without exe (e.g. tests and examples).
AC_ARG_ENABLE(onlylib,
AC_HELP_STRING([--enable-onlylib],[Compile the library alone without tests and examples (default is NO)]),
[case "${enableval}" in 
	yes) onlylib=true ;;
	no) onlylib=false ;;
	*) AC_MSG_ERROR(bad_value ${enableval} for --enable-onlylib) ;;
esac],[onlylib=false])
AM_CONDITIONAL(HAVE_ONLYLIB, test x$onlylib = xtrue)

# Whether to compile support for surrogates
AC_ARG_ENABLE(surrog,
AC_HELP_STRING([--enable-surrog],[Compile the library with support for surrogates (default is YES)]),
[case "${enableval}" in 
	yes) surrog=true ;;
	no) surrog=false ;;
	*) AC_MSG_ERROR(bad_value ${enableval} for --enable-surrog) ;;
esac],[surrog=true])
AM_CONDITIONAL(HAVE_SURROG, test x$surrog = xtrue)
if test x"$surrog" = x"true"; then
   AC_DEFINE(HAVE_SURROG, 1, [surrog])
fi

# Whether debug mode
AC_ARG_ENABLE(debug,
AC_HELP_STRING([--enable-debug],[Enable debug mode (default is NO)]),
[case "${enableval}" in 
	yes) debug=true ;;
	no) debug=false ;;
	*) AC_MSG_ERROR(bad_value ${enableval} for --enable-debug) ;;
esac],[debug=false])
if test x"$debug" = x"true"; then
   AC_DEFINE(HAVE_DEBUG, 1, [debug])
fi
#AM_CONDITIONAL(HAVE_DEBUG, test x$debug = xtrue)

AC_OUTPUT(Makefile src/Makefile tests/Makefile examples/Makefile python/Makefile)

