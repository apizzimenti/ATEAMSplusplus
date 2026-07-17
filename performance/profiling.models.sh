#!/bin/zsh

MODELS=("Bernoulli" "Glauber" "SwendsenWang" "Invasion" "InvadedCluster")
SCALES=(2 4 5 8 11 16)
FIELDS=(2 3)
DIMENSIONS=(2 4)
ITERATIONS=10

HOST=$(hostname -f)
DIR="performance"
mkdir -p ./$DIR/profiling

for MODEL in "${MODELS[@]}"; do
	for SCALE in "${SCALES[@]}"; do
		for FIELD in "${FIELDS[@]}"; do
			for DIMENSION in "${DIMENSIONS[@]}"; do
				# Delete existing ones first.
				PADDED=${(l(3)(0))SCALE}
				PREFIX="$HOST.$MODEL.$PADDED.$DIMENSION.$FIELD.$ITERATIONS"

				echo "## $PREFIX"

				# Record basic statistics.
				perf stat -o $DIR/profiling/$PREFIX.stat ./build/profiling.$MODEL $SCALE $DIMENSION $((DIMENSION/2)) $FIELD $ITERATIONS

				# Find bottlenecks.
				perf record --call-graph fp -o $DIR/profiling/$PREFIX.record ./build/profiling.$MODEL $SCALE $DIMENSION $((DIMENSION/2)) $FIELD $ITERATIONS

				# Make a plot?
				perf script --input=$DIR/profiling/$PREFIX.record \
					| c++filt \
					| gprof2dot --format=perf --strip \
					| dot -Tpng -o $DIR/profiling/$PREFIX.png

				echo "\n"
			done
		done
	done
done
