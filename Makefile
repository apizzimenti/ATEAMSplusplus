
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

install:
	sudo cmake --install build -v

.PHONY: clean build

push: FORCE
	@./push.sh -p
	@./push.sh -m

pull: FORCE
	@./pull.sh -p
	@./pull.sh -m


profiling: FORCE
	@rm -f performance/profiling/*(N)
	./performance/profiling.addition.sh &> ./performance/profiling.addition.log &


timing: FORCE
	@rm -f performance/timing/*(N)
	./performance/timing.addition.sh &> ./performance/timing.addition.log &


test:
	@cd build; ctest

##########
## DOCS ##
##########
docs: FORCE
	@sed -e 's|`\$$|\\f$$|g' -e 's|$$`|\\f$$|g' README.md > README.friendly.md
	@doxygen
	@open file:///Users/apizzimenti/Dropbox/School/GMU/Research/Software/ATEAMS++/docs/index.html

FORCE: ;

