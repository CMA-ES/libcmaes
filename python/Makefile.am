lcmaesdir = $(libdir)

TARGET = lcmaes

if !HAVE_CLANG
  CXXFLAGS += -fopenmp
endif

$(TARGET).so: $(TARGET).o
	$(CXX) -fPIC -shared $(LDFLAGS) $(TARGET).o -L../src/.libs -l$(BOOST_PYTHON_LIB) $(PYTHON_LDFLAGS) -lcmaes -o $(TARGET).so

$(TARGET).o: $(TARGET).cc
	$(CXX) -g -Wall $(CXXFLAGS) -I../ -I$(top_srcdir)/src -I$(EIGEN3_INC) $(PYTHON_CPPFLAGS) -I$(NUMPY_INCLUDE_DIR) -fPIC -c $(srcdir)/$(TARGET).cc

all:lcmaes.so

install: all
	$(INSTALL) lcmaes.so $(DESTDIR)$(libdir)

uninstall-hook:
	rm $(libdir)/lcmaes.so

clean:
	rm *.so *.o
