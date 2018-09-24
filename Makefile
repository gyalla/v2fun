# Required libraries for linker
INC:= -I${GSL_INC} -I$(BOOST_INC)
LDFLAGS:= -L${GSL_LIB} -L$(BOOST_LIB)
LDLIBS  := -lboost_program_options -lgsl -lgslcblas
EXEC    := v2fun
EXECDIR := ../

# Look for libraries if invoking target matches install, check, or coverage
#ifneq (,$(findstring ${MAKECMDGOALS},install-check-coverage))
#    # Find GSL
#    ifdef TACC_GSL_LIB
#        $(info GSL found at ${TACC_GSL_DIR})
#        LDFLAGS:=-L${TACC_GSL_LIB}
#        INC:=-I${TACC_GSL_INC}
#    else
#        ifdef GSL_DIR
#            $(info GSL found at ${GSL_DIR})
#            LDFLAGS:=-L${GSL_DIR}/lib
#            INC:=-I${GSL_DIR}/include
#        else
#            $(info Assuming GSL is in /usr/lib and usr/include)
#            LDFLAGS:=-L/usr/lib
#            INC:=-I/usr/include
#        endif
#    endif
#endif
#
# Export variables so check, install, and coverage can all use libraries
export EXEC
export EXECDIR
export LDFLAGS
export INC
export LDLIBS

info:
	@echo "Available make targets:"
	@echo "  all       : build main program in /src/"
	@echo "  check	    : build and run test unit test suite in /test/unit"
	@echo "  coverage  : build tests with coverage option, run lcov, and generate html in /test/unit/lcov_html"
	@echo "  doc	    : build documentation (doxygen page)" 
all:
	$(MAKE) -C ./src/  



check: 
	$(MAKE)	-C ./test/unit
	$(MAKE) -C ./test/unit check

coverage:
	$(MAKE) -C ./test/unit clobber
	$(MAKE) -C ./test/unit CFLAGS="-O0 -g -Wall --coverage" LDFLAGS="${LDFLAGS} --coverage"
	$(MAKE) -C ./test/unit check 
	./include/lcov/bin/lcov -b ./test/unit --directory ./test/unit --directory ./src/ -no-external -c  -o coverage.info
	./include/lcov/bin/genhtml coverage.info -o lcov_html/
	mv coverage.info test/unit/
	mv lcov_html test/unit/
	@echo
	@echo "-------------------------------------------"
	@echo   Now check test/unit/lcov_html/ for results!
	@echo "-------------------------------------------"
love:
	@echo "not war?"

.PHONY: clobber, clean, doc
clobber: clean 
	-$(MAKE) -C ./src/ clobber
	-$(MAKE) -C ./test/unit/ clobber
	-cd doc/doxygen && rm -rf html && rm -rf latex
	-cd ./test/system/ && rm -rf Sol.png
	-rm -rf $(EXEC)
clean: 
	-$(MAKE) -C ./test/unit/ clean
	-$(MAKE) -C ./src/ clean
	

doc:
	cd doc/doxygen/ && doxygen v2f.dox



