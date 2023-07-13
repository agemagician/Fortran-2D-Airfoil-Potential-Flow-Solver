FC = gfortran
STD := -std=f2008 -Wall -Wextra -pedantic -Wno-compare-reals -ffree-line-length-none
REL := $(STD) -O2
DBG := $(STD) -O0 -g -fcheck=all

# dependency order; only files that exist are used
ORDER := kinds constants naca_geometry grid csr_solver potential_flow postproc io_output
SRC := $(foreach m,$(ORDER),$(wildcard src/$(m).f90))
SUP := $(foreach m,test_support tolerances reference_geometry,$(wildcard tests/$(m).f90))
MAIN := $(wildcard src/main.f90)
TESTS := $(foreach t,csr_poisson naca_geometry grid_diagnostics uniform_flow cylinder joukowski naca_cl io_and_config,$(wildcard tests/test_$(t).f90))

.PHONY: all debug test convergence clean

all:
	@mkdir -p build/release
	$(FC) $(REL) -J build/release $(if $(MAIN),$(SRC) $(MAIN) -o build/release/airfoil,-fsyntax-only $(SRC))

debug:
	@mkdir -p build/debug
	$(FC) $(DBG) -J build/debug $(if $(MAIN),$(SRC) $(MAIN) -o build/debug/airfoil,-fsyntax-only $(SRC))

test:
	@set -e; mkdir -p build/debug; \
	for t in $(TESTS); do n=$$(basename $$t .f90); echo "== $$n"; \
	  $(FC) $(DBG) -J build/debug $(SRC) $(SUP) $$t -o build/debug/$$n; build/debug/$$n; done; \
	if [ -f tests/check_main.sh ]; then $(MAKE) --no-print-directory debug; sh tests/check_main.sh; fi

convergence:
	@mkdir -p build/release
	$(FC) $(REL) -J build/release $(SRC) $(SUP) tests/test_grid_convergence.f90 -o build/release/test_grid_convergence
	build/release/test_grid_convergence

clean:
	rm -rf build
