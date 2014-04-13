all:
	ocamlbuild PhyloCSF.native

test: all
	ocamlbuild testSim.native test.native
	./test.native -verbose

clean:
	rm -f *~
	ocamlbuild -clean

.PHONY: all test clean
