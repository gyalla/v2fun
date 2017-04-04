# FILES
EXEC := v2fun
SRC  := $(wildcard *.cpp)
OBJ  := $(patsubst %.cpp,%.o,$(SRC))

# OPTIONS
CC      := g++ 
CFLAGS  :=-O3 -g -p -Wall
# RULES
$(EXEC): $(OBJ)
	$(LINK.o) -o $@ $^ ${LDFLAGS} $(LDLIBS)
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

