#!/bin/zsh

SCALES=(3 4 6 8 11 16 23 32 45 64)
FIELDS=(2 3 5)
DIMENSIONS=(1 2)
ATTEMPTS=1000
EXEC="persistence"

mkdir -p ./profiles/timing

for SCALE in "${SCALES[@]}"; do
	for DIMENSION in "${DIMENSIONS[@]}"; do
		for FIELD in "${FIELDS[@]}"; do
			PADDED=${(l(2)(0))SCALE}
			PREFIX="$EXEC.$SCALE.$DIMENSION.$FIELD.$ATTEMPTS"
			
			echo $PREFIX
			./build/timing.$EXEC $SCALE $DIMENSION $FIELD $ATTEMPTS >> ./profiles/timing/$PREFIX.data
		done
	done
done
