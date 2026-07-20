#!/bin/zsh

EXECS=("addition")
LENGTHS=(240 3840 9375 61440 219615 983040 3513840 15728640 61509375)
# LENGTHS=(219615 983040 3513840 15728640)
TRIALS=${1:-1000}

HOST=$(hostname -f)
DIR="performance"
mkdir -p ./$DIR/timing

for EXEC in "${EXECS[@]}"; do
	# If we've already done this number of trials, delete the data.
	PREFIX="$HOST.$EXEC.$TRIALS"
	rm ./$DIR/$PREFIX.*

	for LENGTH in "${LENGTHS[@]}"; do

		PADDEDLENGTH=${(r(8)())LENGTH}

		./build/timing.$EXEC $TRIALS $LENGTH >> ./$DIR/timing/$PREFIX.tsv
		echo "completed $PADDEDLENGTH"
	done
done