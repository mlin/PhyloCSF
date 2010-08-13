PKGS = batteries,CamlPaml

all:
	ocamlbuild all.otarget \
		-ocamlc "ocamlfind ocamlc -package ${PKGS} -linkpkg" \
		-ocamlopt "ocamlfind ocamlopt -package ${PKGS} -linkpkg" \
		-ocamldep "ocamlfind ocamldep -package ${PKGS}"

clean:
	rm -f *~
	ocamlbuild -clean
