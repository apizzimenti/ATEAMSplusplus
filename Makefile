
#########
## ENV ##
#########
SHELL := /bin/zsh
DIR := $(dir $(abspath $(firstword $(MAKEFILE_LIST))))

clean:
	rm -rf build

reset: clean
	sudo rm -rf /usr/local/include/ATEAMS++ /usr/local/include/ATEAMS++.h

build:
	cmake -B build
	cmake --build build -v

install: build
	sudo cmake --install build -v

.PHONY: clean build

profile: FORCE
	@screen -dmS profiling ./profiles/profiles.sh

test:
	@cd build; ctest

##########
## DOCS ##
##########
docs: FORCE
	@doxygen

FORCE: ;

