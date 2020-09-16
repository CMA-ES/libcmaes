
SUBDIRS = src
if HAVE_ONLYLIB
else
SUBDIRS += tests examples
endif
if HAVE_PYTHON
SUBDIRS += python
endif

libcmaesincludedir=$(includedir)/libcmaes
libcmaesinclude_HEADERS = include/libcmaes/libcmaes_config.h include/libcmaes/cmaes_export.h

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = libcmaes.pc
