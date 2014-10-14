
FC=gfortran

TESTS=$(patsubst %.f95,%.test,$(wildcard tests/test_*.f95))

test: $(TESTS)
	@rm -f tests/*.out; \
	ok=0; \
	for t in $(TESTS); do \
		b="`basename $$t .test`"; \
		c="tests/$$b.cmp"; \
		log="$$t.out"; \
		echo "--------------------------------------------" > "$$log"; \
		echo "$$b" >> "$$log"; \
		echo "--------------------------------------------" >> "$$log"; \
		echo -n "$$b... "; ./$$t >> "$$log" 2>&1; \
		result=$$?; \
		if test "$$result" = "0"; then \
			if grep -q FAIL "$$log"; then result=1; else result=0; fi; \
		fi; \
		if test -f "$$c"; then \
			if diff -b -u "$$log" "$$c" >> "$$log.tmp"; then true; else result=1; fi; \
			mv -f "$$log.tmp" "$$log"; \
		fi; \
		if test "$$result" = "0"; then echo "OK"; rm -f "$$log"; else echo "FAIL"; ok=1; fi; \
	done; \
	for f in tests/*.out; do test -f "$$f" && cat "$$f"; done; \
	exit "$$ok"

%.o: %.f95
	$(FC) -c -o $@ $^

tests/%.test: tests/%.f95 adjac.o
	$(FC) -o $@ -Itests $^

.PHONY: test
