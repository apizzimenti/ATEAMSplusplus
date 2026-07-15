#!/bin/zsh

MODELS=("Bernoulli" "Glauber" "InvadedCluster" "Invasion" "SwendsenWang")
SCALES=(2 4 5 8 11 16 22 32)
FIELDS=(2 3)
TOPDIMENSIONS=(2 4)
ITERATIONS=1000

DIR="performance"
mkdir -p ./$DIR/profiling

for MODEL in "${MODELS[@]}"; do
	for FIELD in "${FIELDS[@]}"; do
		for SCALE in "${SCALES[@]}"; do
			for TOPDIMENSION in "${TOPDIMENSIONS[@]}"; do
				# Delete existing ones first.
				PADDED=${(l(2)(0))SCALE}
				PREFIX="$MODEL.$FIELD.$PADDED.$TOPDIMENSION"

				echo "########################################"
				echo "## $PREFIX"
				echo "########################################"

				# Record basic statistics.
				perf stat -o $DIR/profiling/$PREFIX.stat ./build/profiling.$MODEL $SCALE $TOPDIMENSION $((TOPDIMENSION/2)) $FIELD $ITERATIONS

				# Find bottlenecks.
				perf record -g -o $DIR/profiling/$PREFIX.record ./build/profiling.$MODEL $SCALE $TOPDIMENSION $((TOPDIMENSION/2)) $FIELD $ITERATIONS

				# Make a plot?
				perf script --input=$DIR/profiling/$PREFIX.record \
					| c++filt \
					| gprof2dot --format=perf --depth=3 --strip --node-thres=5 \
					| dot -Tpng -o $DIR/profiling/$PREFIX.png

				echo "\n"
			done
		done
	done
done
