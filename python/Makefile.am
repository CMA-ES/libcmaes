TARGET = lcmaes

$(TARGET).so: $(TARGET).o
	g++ -fPIC -shared -Wl,--export-dynamic $(TARGET).o -L../src/.libs -l$(BOOST_PYTHON_LIB) $(PYTHON_LDFLAGS) -lcmaes -o $(TARGET).so

$(TARGET).o: $(TARGET).cc
	g++ -g -Wall -std=gnu++11 -fopenmp -I../ -I../src -I$(EIGEN3_INC) $(PYTHON_CPPFLAGS) -fPIC -c $(TARGET).cc

all:lcmaes.so

clean:
	rm *.so *.o
