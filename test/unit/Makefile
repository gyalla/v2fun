# FILES
EXEC := test
SRC  := $(wildcard *.cpp)
COV1 := $(wildcard *.gcno)
COV2 := $(wildcard *.gcda)
OBJ  := $(patsubst %.cpp,%.o,$(SRC))

# OPTIONS
CC      := g++ 
CFLAGS  := -O0 -g -Wall 
OTHER   := ../../src/finiteDiff.cpp \
           ../../src/setup.cpp      \
           ../../src/computeTerms.cpp \
           ../../src/systemSolve.cpp \
           ../../src/Grid.cpp
# RULES


$(EXEC): $(OBJ)
	$(LINK.o) $(OTHER) -o $@ $^ $(INC) $(LDFLAGS) $(LDLIBS)
%.o: %.cpp
	$(COMPILE.c)  $< -o $@ $(INC) $(CFLAGS)

check: 
	$(RM) test_output.txt
	./test

.PHONY: clobber clean
clobber:
	-$(RM) $(EXEC)
	-rm -rf lcov_html
	-rm -f coverage.info
clean: 
	-$(RM) $(OBJ)
	-$(RM) $(COV1)
	-$(RM) $(COV2)
	
