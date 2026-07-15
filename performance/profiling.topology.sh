#!/bin/zsh

EXECS=("persistence")
SCALES=(2 4 5 8 11 16 22 32 45 64)
FIELDS=(3)
DIMENSIONS=(2 4)
ATTEMPTS=4

DIR="performance"
mkdir -p ./$DIR/profiling

for EXEC in "${EXECS[@]}"; do
	for SCALE in "${SCALES[@]}"; do
		for DIMENSION in "${DIMENSIONS[@]}"; do
			for FIELD in "${FIELDS[@]}"; do
				PADDED=${(l(2)(0))SCALE}
				PREFIX="$EXEC.$PADDED.$DIMENSION.$FIELD.$ATTEMPTS"

				echo "## $PREFIX"

				# Record basic statistics.
				perf stat -o $DIR/profiling/$PREFIX.stat ./build/profiling.$EXEC $SCALE $((DIMENSION/2)) $FIELD $ATTEMPTS

				# Find bottlenecks.
				perf record -g -o $DIR/profiling/$PREFIX.record ./build/profiling.$EXEC $SCALE $((DIMENSION/2)) $FIELD $ATTEMPTS

				# Make a plot?
				perf script --input=$DIR/profiling/$PREFIX.record \
					| c++filt \
					| gprof2dot --format=perf --depth=5 --strip --node-thres=5 \
					| dot -Tpng -o $DIR/profiling/$PREFIX.png

				echo "\n"
			done
		done
	done
done
