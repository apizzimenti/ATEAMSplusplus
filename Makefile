
clean:
	rm -rf build

build:
	cmake -B build
	cmake --build build -v

test: clean build
	@cd build; ctest

docs: FORCE
	@doxygen

.PHONY: clean build

FORCE: ;

