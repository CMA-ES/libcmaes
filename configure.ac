# configure for libcmaes
# author
AC_INIT(libcmaes,0.5)
AC_CONFIG_SRCDIR(Makefile.am)
AM_INIT_AUTOMAKE
AM_CONFIG_HEADER(config.h)

ac_default_prefix=`pwd`
bindir=$ac_default_prefix/bin/$host
libdir=$ac_default_prefix/lib/$host

AC_PROG_INSTALL
AC_CANONICAL_HOST

AC_PROG_CC(g++ c++) # force c++ compiler on C files, otherwise BBOB lib doesn't link properly
AC_PROG_CXX(c++ g++ gcc)
AC_PROG_INSTALL
AC_PROG_LIBTOOL

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
SAVE_CFLAGS="$CFLAGS"
SAVE_LIBS="$LIBS"
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

# Check if there is google-glog library installed.
SAVE_CLOG="$CLOG"
SAVE_LIBS="$LIBS"
AC_ARG_WITH(glog, AS_HELP_STRING[--with-glog=GLOG_DIR],
  GLOG_CLOG="-I${with_glog}/include"
  GLOG_LIBS="-L${with_glog}/lib -lglog"
  CLOG="$CLOG $GLOG_CLOG"
  LIBS="$LIBS $GLOG_LIBS"
)
AC_CHECK_LIB(glog, main, ac_cv_have_libglog=1, ac_cv_have_libglog=0)
if test x"$ac_cv_have_libglog" = x"1"; then
   AC_DEFINE(HAVE_LIB_GLOG, 1, [define if you have google glog library])
   if test x"$GLOG_LIBS" = x""; then
      GLOG_LIBS="-lglog"
   fi
else
  GLOG_CLOG=
  GLOG_LIBS=
fi
AC_SUBST(GLOG_LIBS)

# Whether to compile with BBOB
AC_ARG_ENABLE(bbob,
AC_HELP_STRING([--enable-bbob],[Enable the linear programming optimization (default is NO)]),
[case "${enableval}" in 
	yes) bbob=true ;;
	no) bbob=false ;;
	*) AC_MSG_ERROR(bad_value ${enableval} for --enable-bbob) ;;
esac],[bbob=false])
AM_CONDITIONAL(HAVE_BBOB, test x$bbob = xtrue)

AC_OUTPUT(Makefile src/Makefile tests/Makefile examples/Makefile)

