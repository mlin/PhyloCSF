all:
	ocamlbuild all.otarget
	
test:
	ocamlbuild test.otarget
	_build/TestApproxQ.native -verbose

clean:
	rm -f *~
	ocamlbuild -clean
