
clean:
	rm -rf build

build: clean
	cmake -B build
	cmake --build build -v

.PHONY: clean build

