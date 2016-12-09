all: info

info:
	@echo "Available make targets:"
	@echo "  install   : build main program in /src/"
	@echo "  check     : build and run test unit test suite in /test/unit"
	@echo "  coverage  : build tests w/ coverage option, run lcov, and generate html in /test/unit/lcov_html"
	@echo "  doc       : build documentation (doxygen page, and writeup)" 

install:
	$(MAKE) -C ./src/  

check: 
	$(MAKE)	-C ./test/unit
	$(MAKE) -C ./test/unit check

coverage:
	$(MAKE) -C ./test/unit clobber
	$(MAKE) -C ./test/unit CFLAGS="-O0 -g -Wall --coverage" LDFLAGS="-L$$TACC_GSL_LIB -L$$TACC_GRVY_LIB --coverage"
	$(MAKE) -C ./test/unit check 
	./include/lcov/bin/lcov -b ./test/unit --directory ./test/unit --directory ./src/ -no-external -c  -o coverage.info
	./include/lcov/bin/genhtml coverage.info -o lcov_html/
	mv coverage.info test/unit/
	mv lcov_html test/unit/
love:
	@echo "not war?"

.PHONY: clobber, clean, doc
clobber: clean 
	-$(MAKE) -C ./src/ clobber
	-$(MAKE) -C ./test/unit/ clobber
	-$(MAKE) -C ./doc/writeup clean
	-cd doc/doxygen && rm -rf html && rm -rf latex
	-cd ./test/system/ && rm -rf Sol.png

clean: 
	-$(MAKE) -C ./test/unit/ clean
	-$(MAKE) -C ./src/ clean
	-$(MAKE) -C ./doc/writeup clean
	

doc:
	$(MAKE) -C ./doc/writeup/
	cd doc/doxygen/ && doxygen v2f.dox



