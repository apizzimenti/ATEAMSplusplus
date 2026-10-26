#!/bin/zsh

MODELS=("SwendsenWang")				# Models to test.
SCALES=(3 4 5 8 11 16)				# Linear scale of the discrete lattice.
DIMENSIONS=(4)						# Top-level dimension of the discrete lattice.
FIELDS=(3)							# Finite field.
TRIALS=${1:-100}					# Number of trials; default to 100.
SCHEMAS=(1)							# Parallel block allocation schema; default 2*threads.
PROFILINGCONFIGURATION="pmatmul"	# Profiling configuration.

HOST=$(hostname -f)

# Test persistence.
for MODEL in "${MODELS[@]}"; do
	for DIMENSION in "${DIMENSIONS[@]}"; do
		for FIELD in "${FIELDS[@]}"; do
			for SCALE in "${SCALES[@]}"; do
				for SCHEMA in "${SCHEMAS[@]}"; do
					PADDEDSCALE=${(l(2)(0))SCALE}
					PREFIX="$HOST.$MODEL.$PADDEDSCALE.$DIMENSION.$FIELD.$SCHEMA.$PROFILINGCONFIGURATION"

					STAT="./performance/profiling/$PREFIX.stat"
					RECORD="./performance/profiling/$PREFIX.record"
					IMG="./performance/profiling/$PREFIX.png"
					FILES=($STAT $RECORD $IMG)

					for FILE in "${FILES[@]}"; do
						touch $FILE && rm $FILE
					done

					# If the files exist already, delete them.

					perf stat -o ./performance/profiling/$PREFIX.stat ./build/profiling.$MODEL $HOST $SCALE $DIMENSION $FIELD $TRIALS $SCHEMA
					perf record --call-graph fp -s -o ./performance/profiling/$PREFIX.record ./build/profiling.$MODEL $HOST $SCALE $DIMENSION $FIELD $TRIALS $SCHEMA

					perf script --input=./performance/profiling/$PREFIX.record \
						| c++filt \
						| gprof2dot --format=perf --strip \
						| dot -Tpng -o ./performance/profiling/$PREFIX.png > /dev/null

					
					# ./build/timing.$EXEC $HOST $SCALE $DIMENSION $FIELD $TRIALS $SCHEMA
					# gprof ./build/profiling.$EXEC ./gmon.out > ./$DIR/profiling/$PREFIX.gprof
					# gprof ./build/profiling.$EXEC | gprof2dot --strip $ROOT | dot -Tpng -o ./$DIR/profiling/$PREFIX.png

					echo "completed $PREFIX"
					# echo $?
				done
			done
		done
	done
done