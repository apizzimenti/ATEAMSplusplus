#!/bin/zsh

EXECS=("persistence")
SCALES=(11 16 22)
FIELDS=(3)
DIMENSIONS=(4)
PARALLELS=(0 1)
ATTEMPTS=1

HOST=$(hostname -f)
DIR="performance"
mkdir -p ./$DIR/profiling

for PARALLEL in "${PARALLELS[@]}"; do
	if [[ $PARALLEL == 1 ]]; then
		SUFFIX="PARALLEL"
	else
		SUFFIX="SERIAL"
	fi

	for EXEC in "${EXECS[@]}"; do
		for SCALE in "${SCALES[@]}"; do
			for DIMENSION in "${DIMENSIONS[@]}"; do
				for FIELD in "${FIELDS[@]}"; do
					PADDED=${(l(3)(0))SCALE}
					PREFIX="$HOST.$EXEC.$PADDED.$DIMENSION.$FIELD.$ATTEMPTS.$SUFFIX"

					echo "## $PREFIX"

					# Record basic statistics.
					perf stat -o $DIR/profiling/$PREFIX.stat ./build/profiling.$EXEC $SCALE $((DIMENSION/2)) $FIELD $ATTEMPTS $PARALLEL

					# Find bottlenecks.
					perf record --call-graph fp -o $DIR/profiling/$PREFIX.record ./build/profiling.$EXEC $SCALE $((DIMENSION/2)) $FIELD $ATTEMPTS $PARALLEL

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
done