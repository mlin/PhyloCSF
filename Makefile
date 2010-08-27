all: PhyloCSF

.PHONY: PhyloCSF CamlPaml SSE2 clean

ARCH := `uname`.`uname -p`+SSE2
export ARCH

PhyloCSF: CamlPaml
	cd PhyloCSF; $(MAKE) clean; $(MAKE) $(MFLAGS)
	cp PhyloCSF/_build/PhyloCSF.native PhyloCSF.$(ARCH)

CamlPaml: SSE2
	cd CamlPaml; $(MAKE) $(MFLAGS) reinstall

SSE2:
	cd SSE2; $(MAKE) $(MFLAGS) reinstall

clean:
	cd SSE2; $(MAKE) clean
	cd CamlPaml; $(MAKE) clean
	cd PhyloCSF; $(MAKE) clean
	rm -f PhyloCSF.*
