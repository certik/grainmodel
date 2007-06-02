INTERFACES = libmeshpy.i
WRAPPERS   = $(INTERFACES:.i=_wrap.cxx)
WRAPPERS2  = $(INTERFACES:.i=_wrap.h)
PROXIES    = $(INTERFACES:.i=.py      )

all: $(WRAPPERS) libmeshpy.cxx libmeshpy.h
	./setup.py build

%_wrap.cxx: %.i numpy.i libmeshpy.h
	swig -c++ -python $<

clean:
	$(RM) -r build
	$(RM) $(WRAPPERS) $(WRAPPERS2)
	$(RM) $(PROXIES)
