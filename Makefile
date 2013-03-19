all: PhyloCSF

.PHONY: PhyloCSF CamlPaml clean

ARCH := $(shell uname).$(shell uname -m)
PHYLOCSF_BASE := $(shell pwd)
export PHYLOCSF_BASE

PhyloCSF: CamlPaml
	$(MAKE) -C src clean
	$(MAKE) -C src $(MFLAGS)
	cp src/_build/PhyloCSF.native PhyloCSF.$(ARCH)

CamlPaml:
	$(MAKE) -C lib/CamlPaml $(MFLAGS) reinstall
	
clean:
	$(MAKE) -C lib/CamlPaml clean
	$(MAKE) -C src clean
	rm -f PhyloCSF.*
