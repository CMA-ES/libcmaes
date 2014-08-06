# location of the Python header files

PYTHON_VERSION = 2.7
PYTHON_INCLUDE = /usr/include/python$(PYTHON_VERSION)

# location of the Boost Python include files and library

BOOST_INC = /usr/include
BOOST_LIB = /usr/lib

# compile mesh classes
TARGET = lcmaes

$(TARGET).so: $(TARGET).o
	g++ -shared -Wl,--export-dynamic $(TARGET).o -L../src/.libs -L$(BOOST_LIB) -lboost_python -L/usr/lib/python$(PYTHON_VERSION)/config -lpython$(PYTHON_VERSION) -lcmaes -o $(TARGET).so

$(TARGET).o: $(TARGET).cc
	g++ -g -Wall -std=gnu++11 -fopenmp -I../ -I../src -I$(EIGEN3_INC) -I$(PYTHON_INCLUDE) -I$(BOOST_INC) -fPIC -c $(TARGET).cc

all:lcmaes.so

clean:
	rm *.so *.o

#liblcmaesdir = $(libdir)

#liblcmaes_LTLIBRARIES=liblcmaes.la
#liblcmaes_la_SOURCES=lcmaes.cc

#bin_PROGRAMS=lcmaes
#lcmaes_SOURCES=lcmaes.cc

#AM_CXXFLAGS=-I../src/ -I$(EIGEN3_INC) -I/usr/include -I/usr/include/python2.7 -fpic
#AM_LDFLAGS=-shared -Wl,-soname,"lcmaes.so" -L../src
#LDADD=-lpython2.7 -lboost_python -lcmaes -lm $(GLOG_LIBS) $(GFLAGS_LIBS)
