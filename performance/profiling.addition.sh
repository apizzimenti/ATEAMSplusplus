#!/bin/zsh

EXECS=("addition")
LENGTHS=(3513840 15728640)
TRIALS=${1:-100}
PARALLELS=(1)

HOST=$(hostname -f)
DIR="performance"
mkdir -p ./$DIR/timing

PARALLELROOT="ATEAMS::arithmetic::parallelSparseVectorAddition<ATEAMS::Zp>(SparseRREF::sparse_vec<ATEAMS::Zp::dtype, unsigned int>&, SparseRREF::sparse_vec<ATEAMS::Zp::dtype, unsigned int>&, ATEAMS::Ring*, ATEAMS::arithmetic::ComputeOptions<ATEAMS::Zp>&)::{lambda()#1}::operator()() const"
SERIALROOT="ATEAMS::arithmetic::serialSparseVectorAddition<ATEAMS::Zp>(SparseRREF::sparse_vec<ATEAMS::Zp::dtype, unsigned int>&, SparseRREF::sparse_vec<ATEAMS::Zp::dtype, unsigned int>&, ATEAMS::Ring*)"

for EXEC in "${EXECS[@]}"; do
	for LENGTH in "${LENGTHS[@]}"; do
		for PARALLEL in "${PARALLELS[@]}"; do
			# If we've already done this number of trials, delete the data.
			PP="serial"
			ROOT=""
			if [[ $PARALLEL > 0 ]]; then
				PP="parallel"
				ROOT="--root=$PARALLELROOT"
			fi

			PADDEDLENGTH=${(l(8)(0))LENGTH}
			# PADDEDLENGTH=$LENGTH
			PREFIX="$HOST.$EXEC.$TRIALS.$PP.$PADDEDLENGTH"

			# # Record basic statistics.
			# perf stat -o $DIR/profiling/$PREFIX.stat ./build/profiling.$EXEC $HOST $TRIALS $LENGTH $PARALLEL

			# # Find bottlenecks.
			# perf record --call-graph fp -o $DIR/profiling/$PREFIX.record ./build/profiling.$EXEC $HOST $TRIALS $LENGTH $PARALLEL
			./build/profiling.$EXEC $HOST $TRIALS $LENGTH $PARALLEL &&
			gprof ./build/profiling.$EXEC ./gmon.out > ./$DIR/profiling/$PREFIX.gprof &&
			gprof ./build/profiling.$EXEC | gprof2dot --strip $ROOT | dot -Tpng -o ./$DIR/profiling/$PREFIX.gprof.png &&
			
			echo "completed $PADDEDLENGTH"
		done
	done
done