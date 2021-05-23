LIBDIR  = lib
TESTDIR = test

.PHONY: all
all: quantlib test-build

.PHONY: clean
clean:
	cd $(LIBDIR) && $(MAKE) clean
	cd $(TESTDIR) && $(MAKE) clean

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
