CXX=g++
CXXFLAGS=-W -Wall -Wfatal-errors
INCLUDE=-I/usr/include/python2.4

all: py

py:
	swig -c++ -python libmeshpy.i
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c libmeshpy_wrap.cxx
	$(CXX) -shared libmeshpy_wrap.o -o _libmeshpy.so

clean:
	@rm -rf libmeshpy_wrap.cxx libmeshpy_wrap.o libmeshpy.so libmeshpy.py
