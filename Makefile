all: PhyloCSF

.PHONY: PhyloCSF CamlPaml clean

ARCH := $(shell uname).$(shell uname -m)

PhyloCSF: CamlPaml
	$(MAKE) -C src clean
	$(MAKE) -C src $(MFLAGS)
	cp src/_build/PhyloCSF.native PhyloCSF.$(ARCH)

CamlPaml:
	$(MAKE) -C lib/CamlPaml $(MFLAGS) reinstall

test: PhyloCSF
	$(MAKE) -C lib/CamlPaml $(MFLAGS) test
	$(MAKE) -C src $(MFLAGS) test

clean:
	$(MAKE) -C lib/CamlPaml clean
	$(MAKE) -C src clean
	rm -f PhyloCSF.*
