FC=gfortran
CXX=g++

CXXFLAGS=-g -O3
FFLAGS=-g -O3
#FFLAGS=-ggdb -Og

TESTS=$(patsubst %.f95,%.test,$(wildcard tests/test_*.f95)) $(patsubst %.f95,%.test_tapeless,$(wildcard tests/test_*.f95))
EXAMPLES=$(patsubst %.f95,%,$(wildcard examples/*.f95)) $(patsubst %.f95,%_tapeless,$(wildcard examples/*.f95))
CXXEXAMPLES=$(patsubst %.cpp,%,$(wildcard examples/*.cpp))

ADOLC_CFLAGS=-I/usr/include/adolc
ADOLC_LIBS=-ladolc
ADEPT_CFLAGS=-Iadept-1.0/include
ADEPT_LIBS=-Ladept-1.0/lib -ladept
CPPAD_CFLAGS=$(shell pkg-config --cflags cppad)
CPPAD_LIBS=$(shell pkg-config --libs cppad)

all: examples test libadjac.a libadjac_tapeless.a

examples: $(EXAMPLES)

test: $(TESTS)
	@rm -f tests/*.out; \
	ok=0; \
	for t in $(TESTS); do \
		b="`basename $$t .test`"; \
		b="`basename $$b .test_tapeless`"; \
		c="$$t.cmp"; \
		if test ! -f "$$c"; then c="tests/$$b.cmp"; fi; \
		log="$$t.out"; \
		echo "--------------------------------------------" > "$$log"; \
		echo "$$b" >> "$$log"; \
		echo "--------------------------------------------" >> "$$log"; \
		echo -n "$$t... "; ./$$t >> "$$log" 2>&1; \
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
	python generate.py $< $@

adjac_tapeless.f95: adjac.f95.in generate.py
	python generate.py -DTAPELESS=True $< $@

adjac_fft.f95: adjac_fft.f95.in generate.py
	python generate.py $< $@

fftpack/%.f95: fftpack/%.f95.in generate.py
	python generate.py $< $@

sparse_sum.c: sparse_sum.c.in generate.py
	python generate.py sparse_sum.c.in sparse_sum.c

%.o: %.f95
	@install -d build/base
	$(FC) $(FFLAGS) -Jbuild/base -c -o $@ $^

adjac_tapeless.o: adjac_tapeless.f95
	@install -d build/tapeless
	$(FC) $(FFLAGS) -Jbuild/tapeless -c -o $@ $^

%.o: %.c
	gcc -std=c99 $(FFLAGS) -c -o $@ $^

libadjac.a: adjac.o sparse_sum.o build/base/adjac_fft.o build/base/zfftf1.o build/base/zfftb1.o build/base/zffti1.o
	ar cru $@ $^

libadjac_tapeless.a: adjac_tapeless.o sparse_sum.o build/tapeless/adjac_fft.o build/tapeless/zfftf1.o build/tapeless/zfftb1.o build/tapeless/zffti1.o
	ar cru $@ $^

build/base/adjac_fft.o: adjac_fft.f95 adjac.o
	$(FC) $(FFLAGS) -Jbuild/base -c -o $@ $<

build/base/zfftf1.o: fftpack/zfftf1.f95 adjac.o
	$(FC) $(FFLAGS) -Jbuild/base -c -o $@ $<

build/base/zfftb1.o: fftpack/zfftb1.f95 adjac.o
	$(FC) $(FFLAGS) -Jbuild/base -c -o $@ $<

build/base/zffti1.o: fftpack/zffti1.f95
	$(FC) $(FFLAGS) -Jbuild/base -c -o $@ $<

build/tapeless/adjac_fft.o: adjac_fft.f95 adjac_tapeless.o
	$(FC) $(FFLAGS) -Jbuild/tapeless -c -o $@ $<

build/tapeless/zfftf1.o: fftpack/zfftf1.f95 adjac_tapeless.o
	$(FC) $(FFLAGS) -Jbuild/tapeless -c -o $@ $<

build/tapeless/zfftb1.o: fftpack/zfftb1.f95 adjac_tapeless.o
	$(FC) $(FFLAGS) -Jbuild/tapeless -c -o $@ $<

build/tapeless/zffti1.o: fftpack/zffti1.f95
	$(FC) $(FFLAGS) -Jbuild/tapeless -c -o $@ $<

tests/%.test: tests/%.f95 libadjac.a
	$(FC) $(FFLAGS) -Jbuild/base -o $@ -Itests $^ -L. -ladjac

tests/%.test_tapeless: tests/%.f95 libadjac_tapeless.a
	$(FC) $(FFLAGS) -Jbuild/tapeless -o $@ -Itests $^ -L. -ladjac_tapeless

examples/%: examples/%.f95 libadjac.a
	$(FC) $(FFLAGS) -Jbuild/base -o $@ $^ -L. -ladjac

examples/%_tapeless: examples/%.f95 libadjac_tapeless.a
	$(FC) $(FFLAGS) -Jbuild/tapeless -o $@ $^ -L. -ladjac_tapeless

examples/%_adolc: examples/%_adolc.cpp
	$(CXX) $(CXXFLAGS) $(ADOLC_CFLAGS) -o $@ $^ $(ADOLC_LIBS)

examples/bench_simple_tapeless_adolc: examples/bench_simple_tapeless_adolc.cpp examples/bench_simple_adolc.cpp
	if pkg-config --atleast-version=2.5 adolc; then \
	    exec $(CXX) $(CXXFLAGS) $(ADOLC_CFLAGS) -o $@ $< $(ADOLC_LIBS); \
	else \
	    exec $(CXX) $(CXXFLAGS) $(ADOLC_CFLAGS) -DOLD_TAPELESS -o $@ $< $(ADOLC_LIBS); \
	fi

examples/%_adept: examples/%_adept.cpp
	$(CXX) $(CXXFLAGS) $(ADEPT_CFLAGS) -o $@ $^ $(ADEPT_LIBS)

examples/%_cppad: examples/%_cppad.cpp
	$(CXX) $(CXXFLAGS) $(CPPAD_CFLAGS) -o $@ $^ $(CPPAD_LIBS)

compare_adolc: $(EXAMPLES) \
	       examples/bench_simple_adolc examples/bench_simple_tapeless_adolc examples/bench_sparse_adolc
	@echo ""
	@echo "-- bench_simple ----------------------------------------"
	@echo "* ADOLC (tape+eval)"
	time ./examples/bench_simple_adolc
	@echo "* ADOLC (tapeless)"
	time ./examples/bench_simple_tapeless_adolc
	@echo "* ADJAC"
	time ./examples/bench_simple
	@echo "* ADJAC (tapeless)"
	time ./examples/bench_simple_tapeless
	@echo ""
	@echo "-- bench_sparse ----------------------------------------"
	@echo "* ADOLC (tape+eval)"
	time ./examples/bench_sparse_adolc
	@echo "* ADJAC"
	time ./examples/bench_sparse
	@echo "* ADJAC (tapeless)"
	time ./examples/bench_sparse_tapeless

compare_adept: $(EXAMPLES) examples/bench_simple_adept examples/bench_advection_adept
	@echo ""
	@echo "-- bench_simple ----------------------------------------"
	@echo "* ADEPT"
	time ./examples/bench_simple_adept
	@echo "* ADJAC"
	time ./examples/bench_simple
	@echo "* ADJAC (tapeless)"
	time ./examples/bench_simple_tapeless
	@echo ""
	@echo "-- bench_advection ----------------------------------------"
	@echo "* ADEPT"
	time ./examples/bench_advection_adept
	@echo "* ADJAC"
	time ./examples/bench_advection
	@echo "* ADJAC (tapeless)"
	time ./examples/bench_advection_tapeless

compare_cppad: $(EXAMPLES) examples/bench_simple_cppad
	@echo ""
	@echo "-- bench_simple ----------------------------------------"
	@echo "* CPPAD"
	time ./examples/bench_simple_cppad
	@echo "* ADJAC"
	time ./examples/bench_simple
	@echo "* ADJAC (tapeless)"
	time ./examples/bench_simple_tapeless

compare_numdiff: $(EXAMPLES) examples/bench_simple_numdiff
	@echo ""
	@echo "-- bench_simple ----------------------------------------"
	@echo "* Numerical differentiation"
	time ./examples/bench_simple_numdiff
	@echo "* ADJAC"
	time ./examples/bench_simple
	@echo "* ADJAC (tapeless)"
	time ./examples/bench_simple_tapeless

clean:
	rm -rf $(EXAMPLES) $(TESTS) build tests/*.out *.o *.a *.mod $(CXXEXAMPLES)

.PHONY: all test examples compare_adolc compare_adept
