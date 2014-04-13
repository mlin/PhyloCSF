OCAMLBUILDFLAGS=
ifdef FORKWORK
OCAMLBUILDFLAGS=-tag pkg_forkwork
endif

all:
	rm -f ForkMaybe.ml
ifdef FORKWORK
	ln -s ForkYes.ml ForkMaybe.ml
else
	ln -s ForkNo.ml ForkMaybe.ml
endif
	ocamlbuild $(OCAMLBUILDFLAGS) PhyloCSF.native

test: all
	ocamlbuild $(OCAMLBUILDFLAGS) testSim.native test.native
	./test.native -verbose

clean:
	rm -f *~
	ocamlbuild -clean

.PHONY: all test clean
