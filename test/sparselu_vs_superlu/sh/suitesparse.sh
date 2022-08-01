#dgterf

# TEST ALL MATRICES AND OUTPUT IN TAB-DELIMITED
NUM_PROCS=30
SIZES=(1024 2048 4096 8192 16384 32768)
DENSITIES=(0.004)
LIBRARIES=(lass superlu)
TEST_DIR=$HOME/workplan-summer-2021-internship-copy/LASs-sparse-npgetrf-and-dense-tile-pivoting/LASs-sparse-npgetrf/test/slass_v_superlu
OUT_FILE=$TEST_DIR/out/compare-004.txt
RB_DIR=$TEST_DIR/matrices/Synthetic

rm $OUT_FILE
echo -e 'Size\tDensity\tLibrary\tRuntime' >> $OUT_FILE

cd "${RB_DIR}"
for SIZE in "${SIZES[@]}"; do
    for LIBRARY in "${LIBRARIES[@]}"; do
        for DENSITY in "${DENSITIES[@]}"; do

            FILE="${SIZE}-${DENSITY}".rb

            # Size
            echo -e -n "${SIZE}\t">> $OUT_FILE

            # Density
            echo -e -n "$DENSITY\t">> $OUT_FILE

            # Library
            echo -e -n "${LIBRARY}\t">> $OUT_FILE

            # Runtime
            numactl --interleave=all $TEST_DIR/main -p $NUM_PROCS -${LIBRARY} < "${RB_DIR}"/"${FILE}" > /dev/null
            echo >> $OUT_FILE
        done
    done
done

cd ../..