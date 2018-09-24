# FILES
EXEC    := v2fun
EXECDIR := ../
SRC     := $(wildcard *.cpp)
OBJ     := $(patsubst %.cpp,%.o,$(SRC))

# OPTIONS
CC      := g++ 
CFLAGS  :=-O3 -g -p -Wall -fopenmp

# RULES
$(EXECDIR)/$(EXEC): $(OBJ)
	$(LINK.o) -p -fopenmp  -o $@ $^ ${LDFLAGS} $(LDLIBS)
%.o: %.cpp
	$(COMPILE.c) $< -fopenmp -o $@ $(INC)

.PHONY: clean 
clean:
	-$(RM) gmon.out
	-$(RM) $(EXEC)
	-rm $(OBJ)

