
FC=gfortran
CXX=g++

CXXFLAGS=-g -O3
FFLAGS=-g -O3
#FFLAGS=-ggdb -Og

TESTS=$(patsubst %.f95,%.test,$(wildcard tests/test_*.f95))
EXAMPLES=$(patsubst %.f95,%,$(wildcard examples/*.f95))
CXXEXAMPLES=$(patsubst %.cpp,%,$(wildcard examples/*.cpp))

ADOLC_CFLAGS=-I/usr/include/adolc
ADOLC_LIBS=-ladolc
ADEPT_CFLAGS=-Iadept-1.0/include
ADEPT_LIBS=-Ladept-1.0/lib -ladept

all: examples test

examples: $(EXAMPLES)

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
		fi; \
		if test "$$result" = "0"; then echo "OK"; rm -f "$$log"; else echo "FAIL"; ok=1; fi; \
		rm -f "$$log.tmp"; \
	done; \
	for f in tests/*.out; do test -f "$$f" && cat "$$f"; done; \
	exit "$$ok"

adjac.f95: adjac.f95.in generate.py
	python generate.py adjac.f95.in adjac.f95

%.o: %.f95
	$(FC) $(FFLAGS) -c -o $@ $^

tests/%.test: tests/%.f95 adjac.o
	$(FC) $(FFLAGS) -o $@ -Itests $^

examples/%: examples/%.f95 adjac.o
	$(FC) $(FFLAGS) -o $@ $^

examples/%_adolc: examples/%_adolc.cpp
	$(CXX) $(CXXFLAGS) $(ADOLC_CFLAGS) -o $@ $^ $(ADOLC_LIBS)

examples/%_adept: examples/%_adept.cpp
	$(CXX) $(CXXFLAGS) $(ADEPT_CFLAGS) -o $@ $^ $(ADEPT_LIBS)

compare_adolc: examples/bench_simple examples/bench_simple_adolc
	@echo "-- bench_simple ----------------------------------------"
	@echo "* ADOLC (tape+eval)"
	time ./examples/bench_simple_adolc
	@echo "* ADJAC"
	time ./examples/bench_simple

compare_adept: examples/bench_simple examples/bench_simple_adept
	@echo "-- bench_simple ----------------------------------------"
	@echo "* ADEPT"
	time ./examples/bench_simple_adept
	@echo "* ADJAC"
	time ./examples/bench_simple

clean:
	rm -f $(EXAMPLES) $(TESTS) tests/*.out *.o adjac.f95 *.mod $(CXXEXAMPLES)

.PHONY: all test examples compare_adolc compare_adept
