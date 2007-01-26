include ../paths

all: solver convert

solver: solver.o chkopts
	@${CLINKER} -o solver solver.o ${SLEPC_LIB}
	@${RM} solver.o

convert: convert.o chkopts
	@${CLINKER} -o convert convert.o ${SLEPC_LIB}
	@${RM} convert.o

distclean:
	${RM} solver convert

include ${SLEPC_DIR}/bmake/slepc_common
