LIBNAME = CamlPaml

lib:
	ocamlbuild ${LIBNAME}.cma ${LIBNAME}.cmxa

test:
	ocamlbuild test.native
	_build/test.native -verbose


install: lib
	ocamlfind install ${LIBNAME} META _build/${LIBNAME}.cmi _build/${LIBNAME}.cma _build/${LIBNAME}.cmxa _build/${LIBNAME}.a

uninstall:
	ocamlfind remove ${LIBNAME}

reinstall:
	make uninstall
	make install

clean:
	rm -f *~
	ocamlbuild -clean

doc:
	ocamlbuild ${LIBNAME}.docdir/index.html
