lcmaesdir = $(libdir)

TARGET = lcmaes

if !HAVE_CLANG
  CXXFLAGS += -fopenmp
endif

if HAVE_TBB
TBB_INCS = -I$(TBB_INC) -DHAVE_TBB
TBB_LIBS = -L$(TBB_LIB)
TBB_LINKS = -ltbb
endif

$(TARGET).so: $(TARGET).o
	$(CXX) -fPIC -shared $(LDFLAGS) $(TARGET).o -L../src/.libs $(TBB_LIBS) -l$(BOOST_PYTHON_LIB) $(PYTHON_LDFLAGS) -lcmaes $(TBB_LINKS) -o $(TARGET).so

$(TARGET).o: $(TARGET).cc
	$(CXX) -g -Wall $(CXXFLAGS) -I../ -I$(top_srcdir)/src -I$(EIGEN3_INC) $(PYTHON_CPPFLAGS) -I$(NUMPY_INCLUDE_DIR) $(TBB_INCS) -fPIC -c $(srcdir)/$(TARGET).cc

all:lcmaes.so

install: all
	$(INSTALL) lcmaes.so $(DESTDIR)$(libdir)

uninstall-hook:
	rm $(libdir)/lcmaes.so

clean:
	rm *.so *.o
