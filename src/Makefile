# FILES
EXEC := main
SRC  := $(wildcard *.cpp)
OBJ  := $(patsubst %.cpp,%.o,$(SRC))

# OPTIONS
CC      := g++ 
#CFLAGS  := -O0 -g -Wall 
LDFLAGS := -L$$TACC_GRVY_LIB -L$$TACC_GSL_LIB
LDLIBS  := -lgrvy -lgsl -lgslcblas
INC     := -I$$TACC_GRVY_INC -I$$TACC_GSL_INC

# RULES
$(EXEC): $(OBJ)
	$(LINK.o) $(LDLIBS) $(LDFLAGS) -o $@ $^
%.o: %.cpp
	$(COMPILE.c) $< -o $@ $(INC)

main.o systemSolve.o computeTerms.o finiteDiff.o: setup.h
main.o systemSolve.o: computeTerms.h
systemSolve.o: finiteDiff.h


.PHONY: check clobber clean echo 
clobber: clean
	$(RM) $(EXEC)
clean: 
	rm $(OBJ)
echo: 
	@echo $(OBJ)

