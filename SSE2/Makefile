.PHONY : all clean install uninstall reinstall

all:
	ocamlc -c -ccopt -O3 -ccopt -msse2 impl.c
	ocamlc -c SSE2.mli
	ocamlmklib -o SSE2 impl.o SSE2.ml

clean:
	rm -f *.cm* *.cmo *.cmx *.o *.a *.so

install: all
	ocamlfind install SSE2 META SSE2.mli SSE2.cmi SSE2.cma SSE2.cmxa dllSSE2.so SSE2.a libSSE2.a

uninstall:
	ocamlfind remove SSE2

reinstall:
	make uninstall
	make install
