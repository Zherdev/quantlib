LIBDIR   = lib
TESTDIR  = test
CNDOTDIR = cnot
HADDIR   = hadamard_n
BENCHDIR = benchmarking

.PHONY: all
all: quantlib test-build

.PHONY: clean
clean:
	cd $(LIBDIR) && $(MAKE) clean
	cd $(TESTDIR) && $(MAKE) clean
	cd $(BENCHDIR)/$(HADDIR) && $(MAKE) clean
	cd $(BENCHDIR)/$(CNDOTDIR) && $(MAKE) clean

.PHONY: quantlib
quantlib:
	cd $(LIBDIR) && $(MAKE)

.PHONY: test-build
test-build:
	cd $(TESTDIR) && $(MAKE)

.PHONY: test-run-mpiexec
test-run-mpiexec: quantlib test-build
	cd $(TESTDIR) && $(MAKE) run-mpiexec

.PHONY: test-run-mpisubmit
test-run-mpisubmit: quantlib test-build
	cd $(TESTDIR) && $(MAKE) run-mpisubmit

.PHONY: benchmarking-polus
benchmarking-polus: quantlib
	cd $(BENCHDIR) && bash polus_run.sh
