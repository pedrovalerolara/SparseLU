NUM_PROCS=30
SIZES=(1024 2048 3072 4096 5120 6144 7168 8192 9216 10240 11264 12288 13312 14336 15360)
LIBRARIES=(lass superlu)
NUMACTL=(yes no)
DENSITY=0.005
TEST_DIR=$HOME/workplan-summer-2021-internship-copy/LASs-sparse-npgetrf-and-dense-tile-pivoting/LASs-sparse-npgetrf/test/slass_v_superlu
OUT_FILE=$TEST_DIR/out/numactl.txt
FILE=
RB_DIR=$TEST_DIR/matrices/Synthetic

rm $OUT_FILE
echo -e 'Size\tDensity\tLibrary\tNumaCTL\tRuntime' >> $OUT_FILE

cd "${RB_DIR}"
for SIZE in "${SIZES[@]}"; do
    FILE="${SIZE}".rb
    for LIBRARY in "${LIBRARIES[@]}"; do
        for USE_NUMACTL in "${NUMACTL[@]}"; do

            # Size
            echo -e -n "${SIZE}\t">> $OUT_FILE

            # Density
            echo -e -n "$DENSITY\t">> $OUT_FILE

            # Library
            echo -e -n "${LIBRARY}\t">> $OUT_FILE

            # NumaCTL Runtime
            if [ "$USE_NUMACTL" = "yes" ]
            then
                echo -e -n "yes\t">> $OUT_FILE
                numactl --interleave=all $TEST_DIR/main -p $NUM_PROCS -${LIBRARY} < "${RB_DIR}"/"${FILE}" > /dev/null
            else
                echo -e -n "no\t">> $OUT_FILE
                $TEST_DIR/main -p $NUM_PROCS -${LIBRARY} < "${RB_DIR}"/"${FILE}" > /dev/null
            fi

            echo >> $OUT_FILE
        done
    done
done

cd ../..