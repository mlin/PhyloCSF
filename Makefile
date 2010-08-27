all: PhyloCSF

.PHONY: PhyloCSF CamlPaml clean

ARCH := `uname`.`uname -p`
export ARCH

PhyloCSF: CamlPaml
	cd PhyloCSF; $(MAKE) clean; $(MAKE) $(MFLAGS)
	cp PhyloCSF/_build/PhyloCSF.native PhyloCSF.$(ARCH)

CamlPaml:
	cd CamlPaml; $(MAKE) $(MFLAGS) reinstall

clean:
	cd CamlPaml; $(MAKE) clean
	cd PhyloCSF; $(MAKE) clean
	rm -f PhyloCSF.*
