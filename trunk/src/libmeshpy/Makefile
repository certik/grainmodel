INTERFACES = libmesh_mod.i
WRAPPERS   = $(INTERFACES:.i=_wrap.cxx)
PROXIES    = $(INTERFACES:.i=.py      )

all: $(WRAPPERS) libmesh_mod.cxx libmesh_mod.h
	./setup.py build

%_wrap.cxx: %.i numpy.i libmesh_mod.h
	swig -c++ -python $<

clean:
	$(RM) -r build
	$(RM) $(WRAPPERS)
	$(RM) $(PROXIES)
