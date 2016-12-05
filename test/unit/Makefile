# FILES
EXEC := test
SRC  := $(wildcard *.cpp)
OBJ  := $(patsubst %.cpp,%.o,$(SRC))

# OPTIONS
CC      := g++ 
CFLAGS  := -O3
LDFLAGS := -L$$TACC_GSL_LIB -L$$TACC_GRVY_LIB
LDLIBS  :=  -lgsl -lgslcblas -lgrvy
INC     := -I$$TACC_GSL_INC -I$$TACC_GRVY_INC
OTHER   := ../../src/finiteDiff.cpp ../../src/setup.cpp ../../src/computeTerms.cpp ../../src/systemSolve.cpp
# RULES
$(EXEC): $(OBJ)
	$(LINK.o) $(LDLIBS) $(OTHER) $(INC) $(LDFLAGS) -o $@ $^
%.o: %.cpp
	$(COMPILE.c)  $< -o $@ $(INC) 

test_computeTerms.o test_systemSolve.o: ../../src/finiteDiff.h ../../src/computeTerms.h ../../src/systemSolve.h ../../src/setup.h

test.o: test_computeTerms.h test_systemSolve.h

.PHONY: clobber clean echo 
clobber: clean
	$(RM) $(EXEC)
clean: 
	rm $(OBJ)
echo: 
	@echo $(OBJ)

