# FILES
EXEC := v2fun
SRC  := $(wildcard *.cpp)
OBJ  := $(patsubst %.cpp,%.o,$(SRC))

# OPTIONS
CC      := g++ 
LDLIBS  := -lgrvy -lgsl -lgslcblas
INC     := -I$$TACC_GRVY_INC -I$$TACC_GSL_INC
CFLAGS  :=-O3 -g -p -Wall
LDFLAGS := -p -L$$TACC_GSL_LIB -L$$TACC_GRVY_LIB
# RULES
$(EXEC): $(OBJ)
	$(LINK.o) $(LDLIBS) -o $@ $^
%.o: %.cpp
	$(COMPILE.c) $< -o $@ $(INC)

main.o: systemSolve.h
computeTerms.o setup.o: setup.h 
setup.o: computeTerms.h
computeTerms.o: finiteDiff.h

.PHONY: clobber clean 
clobber:
	-$(RM) gmon.out
	-$(RM) $(EXEC)
clean: 
	-rm $(OBJ)

