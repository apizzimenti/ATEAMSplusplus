#!/bin/zsh

HOST=$(hostname -f)


if [[ $HOST == "pangolin" ]]; then
	./performance/timing.topology.sh
else
	./performance/profiling.topology.sh
fi
