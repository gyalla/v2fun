# FILES
EXEC := main
SRC  := $(wildcard *.cc)
OBJ  := $(patsubst %.cc,%.o,$(SRC))

# OPTIONS
CC      := g++ 
CFLAGS  := -O3
LDFLAGS := -L$$TACC_GRVY_LIB -L$$TACC_GSL_LIB
LDLIBS  := -lgrvy -lgsl -lgslcblas
INC     := -I$$TACC_GRVY_INC -I$$TACC_GSL_INC

# RULES
$(EXEC): $(OBJ)
	$(LINK.o) $(LDLIBS) $(INC) $(LDFLAGS) -o $@ $^
%.o: %.cc
	$(COMPILE.c)  $< -o $@ $(INC)

main.o: setup.h
main.o: computeTerms.h

.PHONY: clobber clean echo 
clobber: clean
	$(RM) $(EXEC)
clean: 
	rm $(OBJ)
echo: 
	@echo $(OBJ)

