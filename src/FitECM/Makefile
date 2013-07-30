all:
	ocamlbuild -use-ocamlfind all.otarget
	
test:
	ocamlbuild -use-ocamlfind test.otarget
	_build/TestApproxQ.native -verbose

clean:
	rm -f *~
	ocamlbuild -clean
