#!/bin/zsh

EXECS=("addition")
LENGTHS=(9375 61440 219615 983040 3513840 15728640)
TRIALS=${1:-1000}
PARALLELS=(0 1)

HOST=$(hostname -f)
DIR="performance"
mkdir -p ./$DIR/timing

for EXEC in "${EXECS[@]}"; do
	for LENGTH in "${LENGTHS[@]}"; do
		for PARALLEL in "${PARALLELS[@]}"; do
			PADDEDLENGTH=${(r(8)(0))LENGTH}

			./build/timing.$EXEC $HOST $TRIALS $LENGTH $PARALLEL &&
			echo "completed $PADDEDLENGTH"
		done
	done
done