include ../paths

solver: solver.o chkopts
	@${CLINKER} -o solver solver.o ${SLEPC_LIB}
	@${RM} solver.o

distclean:
	${RM} solver

include ${SLEPC_DIR}/bmake/slepc_common
