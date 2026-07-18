#!/bin/zsh

EXECS=("persistence")
SCALES=(2 4 5 8 11 16)
WIDTHS=(512 1024 2048 4096 8192 16384 32768 65536 131072)
FIELDS=(3)
DIMENSIONS=(2 4)
PARALLELS=(0 1)
ATTEMPTS=4

HOST=$(hostname -f)
DIR="performance"
mkdir -p ./$DIR/profiling


for EXEC in "${EXECS[@]}"; do
	for DIMENSION in "${DIMENSIONS[@]}"; do
		for FIELD in "${FIELDS[@]}"; do
			for PARALLEL in "${PARALLELS[@]}"; do
				if [[ $PARALLEL == 1 ]]; then
					SUFFIX="PARALLEL"

					for WIDTH in "${WIDTHS[@]}"; do
						for SCALE in "${SCALES[@]}"; do
							PADDEDSCALE=${(l(3)(0))SCALE}
							PADDEDWIDTH=${(l(6)(0))WIDTH}
							PREFIX="$HOST.$EXEC.$FIELD.$DIMENSION.$PADDEDSCALE.$SUFFIX.$PADDEDWIDTH"

							echo "## $PREFIX"

							# Record basic statistics.
							perf stat -o $DIR/profiling/$PREFIX.stat ./build/profiling.$EXEC $SCALE $((DIMENSION/2)) $FIELD $ATTEMPTS $PARALLEL $WIDTH

							# Find bottlenecks.
							perf record --call-graph fp -o $DIR/profiling/$PREFIX.record ./build/profiling.$EXEC $SCALE $((DIMENSION/2)) $FIELD $ATTEMPTS $PARALLEL $WIDTH

							# Make a plot?
							perf script --input=$DIR/profiling/$PREFIX.record \
								| c++filt \
								| gprof2dot --format=perf --strip \
								| dot -Tpng -o $DIR/profiling/$PREFIX.png

							echo "\n"
						done
					done
				else
					SUFFIX="SERIAL"
					for SCALE in "${SCALES[@]}"; do
						PADDED=${(l(3)(0))SCALE}
						PREFIX="$HOST.$EXEC.$FIELD.$DIMENSION.$PADDED.$SUFFIX"

						echo "## $PREFIX"

						# Record basic statistics.
						perf stat -o $DIR/profiling/$PREFIX.stat ./build/profiling.$EXEC $SCALE $((DIMENSION/2)) $FIELD $ATTEMPTS 0 1

						# Find bottlenecks.
						perf record --call-graph fp -o $DIR/profiling/$PREFIX.record ./build/profiling.$EXEC $SCALE $((DIMENSION/2)) $FIELD $ATTEMPTS 0 1

						# Make a plot?
						perf script --input=$DIR/profiling/$PREFIX.record \
							| c++filt \
							| gprof2dot --format=perf --strip \
							| dot -Tpng -o $DIR/profiling/$PREFIX.png

						echo "\n"
					done
				fi
			done
		done
	done
done