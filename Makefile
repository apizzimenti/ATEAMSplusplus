
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

profile: FORCE
	@screen -dmS profiling ./profiles/profiles.models.sh

timing: FORCE
	@screen -dmS timing ./profiles/timing.persistence.sh

test:
	@cd build; ctest

push: FORCE
	@./push.sh -p

pull: FORCE
	@./retrieve.sh -p

##########
## DOCS ##
##########
docs: FORCE
	@sed -e 's|`\$$|\\f$$|g' -e 's|$$`|\\f$$|g' README.md > README.friendly.md
	@doxygen
	@open file:///Users/apizzimenti/Dropbox/School/GMU/Research/Software/ATEAMS++/docs/index.html

FORCE: ;

