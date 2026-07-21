#!/bin/zsh

EXECS=("addition")
LENGTHS=(219615 983040)
TRIALS=${1:-10}
PARALLELS=(1)

HOST=$(hostname -f)
DIR="performance"
mkdir -p ./$DIR/timing

for EXEC in "${EXECS[@]}"; do
	for PARALLEL in "${PARALLELS[@]}"; do
		# If we've already done this number of trials, delete the data.
		PP="serial"
		if [[ $PARALLEL > 0 ]]; then
			PP="parallel"
		fi

		for LENGTH in "${LENGTHS[@]}"; do
			# PADDEDLENGTH=${(l(8)(0))LENGTH}
			PADDEDLENGTH=$LENGTH
			PREFIX="$HOST.$EXEC.$TRIALS.$PP.$PADDEDLENGTH"

			# Record basic statistics.
			perf stat -o $DIR/profiling/$PREFIX.stat ./build/profiling.$EXEC $HOST $TRIALS $LENGTH $PARALLEL

			# Find bottlenecks.
			perf record --call-graph fp -o $DIR/profiling/$PREFIX.record ./build/profiling.$EXEC $HOST $TRIALS $LENGTH $PARALLEL

			# Make a plot?
			perf script --input=$DIR/profiling/$PREFIX.record \
				| c++filt \
				| gprof2dot --format=perf --strip \
				| dot -Tpng -o $DIR/profiling/$PREFIX.png
			
			echo "completed $PADDEDLENGTH"
		done
	done
done