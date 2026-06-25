
#########
## ENV ##
#########
SHELL := /bin/zsh
DIR := $(dir $(abspath $(firstword $(MAKEFILE_LIST))))

clean:
	rm -rf build

build:
	cmake -B build
	cmake --build build -v

install: build
	sudo cmake --install build -v

.PHONY: clean build

test: clean build install
	@cd build; ctest

##########
## DOCS ##
##########
docs: FORCE
	@doxygen

FORCE: ;

