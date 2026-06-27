#!/bin/zsh

MODELS=("Bernoulli" "Glauber" "InvadedCluster" "InvasionPercolation" "SwendsenWang")
SCALES=(3)
FIELDS=(2)
TOPDIMENSIONS=(2)
ITERATIONS=100

EXECPREFIX="../build"

for MODEL in "${MODELS[@]}"; do
	for FIELD in "${FIELDS[@]}"; do
		for SCALE in "${SCALES[@]}"; do
			for TOPDIMENSION in "${TOPDIMENSIONS[@]}"; do
				# Record basic statistics.
				perf stat -o profiles/reports/$MODEL.$FIELD.$SCALE.$TOPDIMENSION.txt ./build/profiles.$MODEL $SCALE $TOPDIMENSION $((TOPDIMENSION/2)) $FIELD $ITERATIONS &
				wait $!

				# Find bottlenecks.
				perf record -o profiles/reports/$MODEL.$FIELD.$SCALE.$TOPDIMENSION.txt ./build/profiles.$MODEL $SCALE $TOPDIMENSION $((TOPDIMENSION/2)) $FIELD $ITERATIONS &
				wait $!
			done
		done
	done
done
