#!/bin/zsh

EXECS=("persistence")
SCALES=(2 4 5 8 11 16 22 32 45 64 90 128 181 256 362 512)
FIELDS=(2 3)
DIMENSIONS=(2 4)
ATTEMPTS=51
PARALLELS=(0 1)

HOST=$(hostname -f)
DIR="performance"
mkdir -p ./$DIR/timing


for EXEC in "${EXECS[@]}"; do
	for DIMENSION in "${DIMENSIONS[@]}"; do
		for FIELD in "${FIELDS[@]}"; do
			for PARALLEL in "${PARALLELS[@]}"; do
				if [[ $PARALLEL == 1 ]]; then
					SUFFIX="PARALLEL"
				else
					SUFFIX="SERIAL"
				fi

				for SCALE in "${SCALES[@]}"; do
					PADDED=${(l(3)(0))SCALE}
					PREFIX="$HOST.$EXEC.$FIELD.$DIMENSION.$PADDED.$SUFFIX"

					echo $PREFIX
					./build/timing.$EXEC $SCALE $((DIMENSION/2)) $FIELD $ATTEMPTS $PARALLEL >> ./$DIR/timing/$PREFIX.data
				done
			done
		done
	done
done