INTERFACES = libmeshpy.i
WRAPPERS   = $(INTERFACES:.i=_wrap.cxx)
PROXIES    = $(INTERFACES:.i=.py      )

all: $(WRAPPERS) libmeshpy.cxx libmeshpy.h
	./setup.py build

%_wrap.cxx: %.i numpy.i libmeshpy.h
	swig -c++ -python $<

clean:
	$(RM) -r build
	$(RM) $(WRAPPERS)
	$(RM) $(PROXIES)
	py.cleanup
