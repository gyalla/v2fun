# FILES
EXEC := main
SRC  := $(wildcard *.cpp)
OBJ  := $(patsubst %.cpp,%.o,$(SRC))

# OPTIONS
CC      := g++ 
CFLAGS  := -O3
LDFLAGS := -L$$TACC_GRVY_LIB -L$$TACC_GSL_LIB
LDLIBS  := -lgrvy -lgsl -lgslcblas
INC     := -I$$TACC_GRVY_INC -I$$TACC_GSL_INC

# RULES
$(EXEC): $(OBJ)
	$(LINK.o) -g -Wall $(LDLIBS) $(INC) $(LDFLAGS) -o $@ $^
%.o: %.cpp
	$(COMPILE.c) -g -Wall  $< -o $@ $(INC)

main.o systemSolve.o computeTerms.o: setup.h
main.o systemSolve.o: computeTerms.h
systemSolve.o: finiteDiff.h


.PHONY: check clobber clean echo 
clobber: clean
	$(RM) $(EXEC)
clean: 
	rm $(OBJ)
echo: 
	@echo $(OBJ)
check:
	../test/unit/test

