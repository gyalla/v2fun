main:
	make -C ./src/

check: 
	make -C ./test/unit
	./test/unit/test

converge:


.PHONY: clobber 
clobber:
	$(RM) /src/main
	$(RM) /test/unit/test


