all: ../resources/dxPhyloCSF

../resources/dxPhyloCSF: _tags dxPhyloCSF.ml
	ocamlbuild -use-ocamlfind dxPhyloCSF.native
	mkdir -p ../resources
	cp -f dxPhyloCSF.native ../resources/dxPhyloCSF

clean:
	ocamlbuild -clean
	rm -f ../resources/dxPhyloCSF
