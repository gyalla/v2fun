# FILES
EXEC := test
SRC  := $(wildcard *.cpp)
COV1 := $(wildcard *.gcno)
COV2 := $(wildcard *.gcda)
OBJ  := $(patsubst %.cpp,%.o,$(SRC))

# OPTIONS
CC      := g++ 
LDFLAGS := -L$$TACC_GSL_LIB -L$$TACC_GRVY_LIB
CFLAGS  := -O0 -g -Wall 
LDLIBS  :=  -lgsl -lgslcblas -lgrvy
INC     := -I$$TACC_GSL_INC -I$$TACC_GRVY_INC
OTHER   := ../../src/finiteDiff.cpp ../../src/setup.cpp ../../src/computeTerms.cpp ../../src/systemSolve.cpp
# RULES
$(EXEC): $(OBJ)
	$(LINK.o) $(LDLIBS) $(OTHER) $(INC) $(LDFLAGS) -o $@ $^
%.o: %.cpp
	$(COMPILE.c)  $< -o $@ $(INC) 

test_computeTerms.o test_systemSolve.o test.o: ../../src/finiteDiff.h ../../src/computeTerms.h ../../src/systemSolve.h ../../src/setup.h ../../src/main.cpp

test.o: test_computeTerms.h test_systemSolve.h 

check: 
	$(RM) test_output.txt
	./test

.PHONY: clobber clean echo 
clobber: clean
	$(RM) $(EXEC)
	rm -rf lcov_html
	rm -rf coverage.info
clean: 
	$(RM) $(OBJ)
	$(RM) $(COV1)
	$(RM) $(COV2)
	
echo: 
	@echo $(OBJ)


