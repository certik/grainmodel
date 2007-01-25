# $Id: Makefile,v 1.12 2003/11/06 20:30:20 benkirk Exp $

include ../paths

target 	   := ./quantum


all:: $(target)

# The location of the mesh library
#meshdir := ~/libmesh/test3/libmesh
#meshdir := ~/libmesh/test4/libmesh
#meshdir := ~/libmesh/test5/libmesh

# include the library options determined by configure.  This will
# set the variables INCLUDE and LIBS that we will need to build and
# link with the library.
include $(meshdir)/Make.common

# include the dependency list
include .depend
main.$(obj-suffix): configfile.h


###############################################################################
# File management.  This is where the source, header, and object files are
# defined

#
# source files
srcfiles 	:= $(wildcard *.C)

#
# object files
objects		:= $(patsubst %.C, %.$(obj-suffix), $(srcfiles))
###############################################################################



.PHONY: clean clobber distclean

###############################################################################
# Target:
#

# Production rules:  how to make the target - depends on library configuration
$(target): $(objects)
	@echo "Linking "$@"..."
	@$(CXX) $(CXXFLAGS) $(objects) -o $@ $(LIBS) $(LDFLAGS)


# Useful rules.
clean:
	@rm -f $(objects) *~
	@rm -f quantum
	@rm -f .depend

clobber:
	@$(MAKE) clean
	@rm -f $(target) out.gmv

distclean:
	@$(MAKE) clobber
	@rm -f *.o *.g.o *.pg.o

run: $(target)
	@echo "***************************************************************"
	@echo "* Running Example " $(target)
	@echo "***************************************************************"
	@echo " "
	@$(target)
	@echo " "
	@echo "***************************************************************"
	@echo "* Done Running Example " $(target)
	@echo "***************************************************************"
#
# Dependencies
#
.depend:
	@$(perl) $(meshdir)/contrib/bin/make_dependencies.pl $(foreach i, $(wildcard $(meshdir)/include/*), -I$(i)) "-S\$$(obj-suffix)" $(srcfiles) > .depend
	@echo "Updated .depend"

###############################################################################
