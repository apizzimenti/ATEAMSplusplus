
clean:
	rm -rf build

build:
	cmake -B build
	cmake --build build -v

install: build
	sudo cmake --install build -v

test: clean build
	@cd build; ctest

docs: FORCE
	@doxygen

.PHONY: clean build

FORCE: ;

