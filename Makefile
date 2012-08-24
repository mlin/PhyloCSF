all: PhyloCSF

.PHONY: PhyloCSF CamlPaml twt cde-package CDE clean

ARCH := $(shell uname).$(shell uname -m)
PHYLOCSF_BASE := $(shell pwd)
export PHYLOCSF_BASE

PhyloCSF: CamlPaml
	cd src; $(MAKE) clean; PATH=$(PATH):$(CURDIR)/twt $(MAKE) $(MFLAGS)
	cp src/_build/PhyloCSF.native PhyloCSF.$(ARCH)

CamlPaml: twt
	cd lib/CamlPaml; PATH=$(PATH):$(CURDIR)/twt $(MAKE) $(MFLAGS) reinstall

twt:
	git submodule update --init
	cd twt && $(MAKE)

cde-package: PhyloCSF CDE
	rm -rf cde-package cde-package.$(ARCH)
	CDE/cde ./PhyloCSF.$(ARCH) 12flies PhyloCSF_Examples/tal-AA.fa 
	CDE/cde ./PhyloCSF.$(ARCH) 29mammals PhyloCSF_Examples/ALDH2.exon5.fa --frames=3 --allScores
	CDE/cde ./PhyloCSF.$(ARCH) 29mammals PhyloCSF_Examples/Aldh2.mRNA.fa --strategy=fixed --frames=3 --orf=ATGStop --minCodons=400 --allScores --removeRefGaps --aa
	mv cde-package cde-package.$(ARCH)
	tar -cf cde-package.$(ARCH).tar cde-package.$(ARCH)

CDE:
	git submodule update --init
	cd CDE && $(MAKE)

clean:
	cd lib/CamlPaml; $(MAKE) clean
	cd src; $(MAKE) clean
	rm -f PhyloCSF.*
	rm -rf cde-package*
	cd twt && $(MAKE) clean || true
	cd CDE && $(MAKE) clean || true
