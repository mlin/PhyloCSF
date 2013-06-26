all: default

default:
	rm -rf resources/PhyloCSF*
	mkdir -p resources
	cp -r ../PhyloCSF ../PhyloCSF.Linux.x86_64 ../PhyloCSF_Parameters resources
	$(MAKE) -C src all

clean:
	$(MAKE) -C src clean

test: all
	test/test.sh

.PHONY: all default clean test
