#!/bin/zsh

EXECS=("addition")
LENGTHS=(100000 1000000 10000000 100000000)
STRATEGIES=(0 1 2 3 4 5)
TRIALS=${1:-100}
PARALLELS=(0 1)

HOST=$(hostname -f)
DIR="performance"
mkdir -p ./$DIR/timing

# Test pa
for EXEC in "${EXECS[@]}"; do
	for LENGTH in "${LENGTHS[@]}"; do
		# Get serial stuff.
		./build/timing.$EXEC $HOST $TRIALS $LENGTH 0 0
		echo "completed $LENGTH (serial)"

		# Get parallel stuff, trying with different strategies.
		for STRATEGY in "${STRATEGIES[@]}"; do
			./build/timing.$EXEC $HOST $TRIALS $LENGTH 1 $STRATEGY
			echo "completed $LENGTH (parallel, $STRATEGY)"
		done
	done
done