PKGS = gsl

lib:
	ocamlbuild -cflags -for-pack,CamlPaml CamlPaml.cma CamlPaml.cmxa \
		-ocamlc "ocamlfind ocamlc -package ${PKGS}" \
		-ocamlopt "ocamlfind ocamlopt -package ${PKGS}" \
		-ocamldep "ocamlfind ocamldep -package ${PKGS}"

TESTPKGS = ${PKGS},oUnit
test:
	ocamlbuild test.native \
		-ocamlc "ocamlfind ocamlc -package ${TESTPKGS} -linkpkg" \
		-ocamlopt "ocamlfind ocamlopt -package ${TESTPKGS} -linkpkg" \
		-ocamldep "ocamlfind ocamldep -package ${TESTPKGS}"
	_build/test.native -verbose


install: lib
	ocamlfind install CamlPaml META _build/CamlPaml.cmi _build/CamlPaml.cma _build/CamlPaml.cmxa _build/CamlPaml.a

uninstall:
	ocamlfind remove CamlPaml

reinstall:
	make uninstall
	make install

clean:
	rm -f *~
	ocamlbuild -clean

doc:
	ocamlbuild CamlPaml.docdir/index.html
