all: PhyloCSF

.PHONY: PhyloCSF CamlPaml clean

ARCH := `uname`.`uname -p`
export ARCH

PhyloCSF: CamlPaml
	cd src; $(MAKE) clean; $(MAKE) $(MFLAGS)
	cp src/_build/PhyloCSF.native PhyloCSF.$(ARCH)

CamlPaml: 
	cd lib/CamlPaml; $(MAKE) $(MFLAGS) reinstall

clean:
	cd lib/CamlPaml; $(MAKE) clean
	cd src; $(MAKE) clean
	rm -f PhyloCSF.*
