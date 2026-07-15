#!/bin/zsh

EXECS=("persistence")
SCALES=(2 4 5 8 11 16 22 32)
FIELDS=(2 3)
DIMENSIONS=(2 4)
ATTEMPTS=100

DIR="performance"
mkdir -p ./$DIR/timing

for EXEC in "${EXECS[@]}"; do
	for SCALE in "${SCALES[@]}"; do
		for DIMENSION in "${DIMENSIONS[@]}"; do
			for FIELD in "${FIELDS[@]}"; do
				PADDED=${(l(2)(0))SCALE}
				PREFIX="$EXEC.$PADDED.$DIMENSION.$FIELD.$ATTEMPTS"

				echo $PREFIX
				./build/timing.$EXEC $SCALE $((DIMENSION/2)) $FIELD $ATTEMPTS >> ./$DIR/timing/$PREFIX.data
			done
		done
	done
done
