#Uncomment the follow two lines and add the root path to GSL and BOOST!
#GSL_DIR   = 
#BOOST_DIR = 

INC:= -I${GSL_DIR}/include/ -I$(BOOST_DIR)/include/
LDFLAGS:= -L${GSL_DIR}/lib/ -L$(BOOST_DIR)/lib/
LDLIBS  := -lboost_program_options -lgsl -lgslcblas
EXEC    := v2fun
EXECDIR := ../

# Export variables 
export EXEC
export EXECDIR
export LDFLAGS
export INC
export LDLIBS

info:
	@echo "Available make targets:"
	@echo "  all       : build main program"
	@echo "  check	    : build and run test unit test suite in /test/unit"
	@echo "  coverage  : build tests with coverage option, run lcov, and generate html in /test/unit/lcov_html"
	@echo "  doc	    : build documentation (doxygen page)" 
	@echo 
	@echo Remember to edit the GSL_DIR and BOOST_DIR variables at the top of Makefile before building! 

all:
	$(MAKE) -C ./src/  


check: 
	$(MAKE)	-C ./test/unit
	$(MAKE) -C ./test/unit check

coverage:
	@echo "-------------------------------------------------------"
	@echo   Note: Must have lcov installed to use coverage feature
	@echo "-------------------------------------------------------"
	$(MAKE) -C ./test/unit clean
	$(MAKE) -C ./test/unit CFLAGS="-O0 -g -Wall --coverage" LDFLAGS="${LDFLAGS} --coverage"
	$(MAKE) -C ./test/unit check 
	lcov -b ./test/unit --directory ./test/unit --directory ./src/ -no-external -c  -o coverage.info
	genhtml coverage.info -o lcov_html/
	mv coverage.info test/unit/
	mv lcov_html test/unit/
	@echo
	@echo "-------------------------------------------"
	@echo   Now check test/unit/lcov_html/ for results!
	@echo "-------------------------------------------"
love:
	@echo "not war?"

.PHONY: clean, doc
clean: 
	-cd doc/doxygen && rm -rf html && rm -rf latex
	-rm -rf $(EXEC)
	-$(MAKE) -C ./test/unit/ clean
	-$(MAKE) -C ./src/ clean
	
doc:
	cd doc/doxygen/ && doxygen v2f.dox
