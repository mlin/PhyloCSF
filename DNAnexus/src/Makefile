all: build install

build:
	mkdir -p build && cd build && cmake -D CMAKE_BUILD_TYPE:STRING=RELEASE .. && $(MAKE)

install: build
	cp -f build/PhyloCSF ../resources

clean:
	rm -rf build/
	rm -f ../resources/PhyloCSF

.PHONY: all clean build install
