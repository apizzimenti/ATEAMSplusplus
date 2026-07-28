
#########
## ENV ##
#########
SHELL := /bin/zsh
DIR := $(dir $(abspath $(firstword $(MAKEFILE_LIST))))
NPROCS := 1
OS := $(shell uname -s)

ifeq ($(OS),Linux)
	NPROCS := $(shell grep -c ^processor /proc/cpuinfo)
endif

ifeq ($(OS),Darwin)
	NPROCS := $(shell sysctl -n hw.ncpu)
endif

clean:
	rm -rf build

reset: clean
	sudo rm -rf /usr/local/include/ATEAMS++ /usr/local/include/ATEAMS++.h

build:
	cmake -B build
	cmake --build build -v -- -j $(NPROCS)
install:
	sudo cmake --install build -v

.PHONY: clean build

push: FORCE
	@./push.sh -p
	@./push.sh -m

pull: FORCE
	@./pull.sh -p
	@./pull.sh -m


samples: FORCE
# 	@rm -f performance/timing/samples.*(N)
	./performance/performance.samples.sh


profiling: FORCE
	screen -dmS profiling.persistence ./performance/profiling.persistence.sh


timing: FORCE
	screen -dmS timing.persistence ./performance/timing.persistence.sh


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

