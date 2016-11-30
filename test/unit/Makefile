# FILES
EXEC := test
SRC  := $(wildcard *.cc)
OBJ  := $(patsubst %.cc,%.o,$(SRC))

# OPTIONS
CC      := g++ 
CFLAGS  := -O3
LDFLAGS := -L$$TACC_GSL_LIB -L$$TACC_GRVY_LIB
LDLIBS  :=  -lgsl -lgslcblas -lgrvy
INC     := -I$$TACC_GSL_INC -I$$TACC_GRVY_INC
OTHER   := ../../src/setup.cc ../../src/computeTerms.cc ../../src/systemSolve.cc
# RULES
$(EXEC): $(OBJ)
	$(LINK.o) $(LDLIBS) $(OTHER) $(INC) $(LDFLAGS) -o $@ $^
%.o: %.cc
	$(COMPILE.c)  $< -o $@ $(INC)

test.o: ../../src/setup.h

.PHONY: clobber clean echo 
clobber: clean
	$(RM) $(EXEC)
clean: 
	rm $(OBJ)
echo: 
	@echo $(OBJ)

