OCAMLBUILDFLAGS=
ifdef FORKWORK
OCAMLBUILDFLAGS=-tag "package(forkwork)"
endif

all:
	rm -f ForkMaybe.ml
ifdef FORKWORK
	ln -s ForkYes.ml ForkMaybe.ml
else
	ln -s ForkNo.ml ForkMaybe.ml
endif
	ocamlbuild -use-ocamlfind $(OCAMLBUILDFLAGS) PhyloCSF.native

test: testexe
	./test.native -verbose

testexe: all
	ocamlbuild -use-ocamlfind $(OCAMLBUILDFLAGS) testSim.native test.native

clean:
	rm -f *~
	ocamlbuild -clean

.PHONY: all test testexe clean
