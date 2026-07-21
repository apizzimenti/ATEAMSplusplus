
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
# 	@screen -dmS profiling.models ./performance/profiling.models.sh
	@screen -dmS profiling.topology ./performance/profiling.topology.sh

timing: FORCE
	./performance/timing.addition.sh &> ./performance/timing.addition.log &

gauntlet: FORCE
# 	rm -f ./performance/profiling/*(N) ./performance/timing/*(N)
# 	screen -dmS gauntlet.topology ./performance/gauntlet.topology.sh
# 	screen -dmS timing.addition.sample ./performance/timing.addition.sh 50
# 	screen -dmS timing.addition.default ./performance/timing.addition.sh
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

