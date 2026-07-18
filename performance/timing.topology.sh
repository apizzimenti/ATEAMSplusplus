#!/bin/zsh

EXECS=("persistence")
SCALES=(2 4 5 8 11 16 22 32 45 64 90 128 181 256 362 512)
WIDTHS=(512 1024 2048 4096 8192 16384 32768 65536 131072)
FIELDS=(3 2)
DIMENSIONS=(2 4)
ATTEMPTS=51
PARALLELS=(1 0)

HOST=$(hostname -f)
DIR="performance"
mkdir -p ./$DIR/timing


for EXEC in "${EXECS[@]}"; do
	for DIMENSION in "${DIMENSIONS[@]}"; do
		for FIELD in "${FIELDS[@]}"; do
			for PARALLEL in "${PARALLELS[@]}"; do
				if [[ $PARALLEL == 1 && $FIELD == 3 ]]; then
					SUFFIX="PARALLEL"

					for SCALE in "${SCALES[@]}"; do
						for WIDTH in "${WIDTHS[@]}"; do
							PADDEDSCALE=${(l(3)(0))SCALE}
							PADDEDWIDTH=${(l(6)(0))WIDTH}
							PREFIX="$HOST.$EXEC.$FIELD.$DIMENSION.$PADDEDSCALE.$SUFFIX.$PADDEDWIDTH"

							echo $PREFIX
							./build/timing.$EXEC $SCALE $((DIMENSION/2)) $FIELD $ATTEMPTS $PARALLEL $WIDTH >> ./$DIR/timing/$PREFIX.data
						done
					done

				else
					SUFFIX="SERIAL"
					for SCALE in "${SCALES[@]}"; do
						PADDED=${(l(3)(0))SCALE}
						PREFIX="$HOST.$EXEC.$FIELD.$DIMENSION.$PADDED.$SUFFIX"

						echo $PREFIX
						./build/timing.$EXEC $SCALE $((DIMENSION/2)) $FIELD $ATTEMPTS 0 1 >> ./$DIR/timing/$PREFIX.data
					done
				fi
			done
		done
	done
done