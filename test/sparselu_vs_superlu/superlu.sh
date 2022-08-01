#dgterf

# TEST ALL MATRICES AND OUTPUT IN TAB-DELIMITED
NUM_PROCS=32
CORES=(25 26 27 28 29 30 31)
SIZES=(1024 2048 3072 4096 5120 6144 7168 8192 9216 10240 11264 12288 13312 14336 15360)
DENSITY=0.005
TEST_DIR=$HOME/workplan-summer-2021-internship-copy/LASs-sparse-npgetrf-and-dense-tile-pivoting/LASs-sparse-npgetrf/test/slass_v_superlu
OUT_FILE=$TEST_DIR/superlu_results.txt
FILE=
RB_DIR=$TEST_DIR/matrices/Synthetic
INDEX=0

rm $OUT_FILE
echo -e 'INDEX\tFILE\tTYPE\tSIZE\tDENSITY\tTHREADS\tRUNTIME_SUPERLU' >> $OUT_FILE

cd "${RB_DIR}"
for SIZE in "${SIZES[@]}"; do
    FILE="${SIZE}".rb
    for NUM_PROCS in "${CORES[@]}"; do
        ((INDEX=INDEX+1))
        echo -e -n $INDEX'\t'>> $OUT_FILE
        echo -e -n $FILE'\t'>> $OUT_FILE
        $TEST_DIR/matrix_attributes "${RB_DIR}"/"${FILE}"
        echo -e -n '\t'>> $OUT_FILE
        echo -e -n $DENSITY>> $OUT_FILE
        echo -e -n '\t'>> $OUT_FILE
        echo -e -n $NUM_PROCS>> $OUT_FILE
        echo -e -n '\t'>> $OUT_FILE
        numactl --interleave=all $TEST_DIR/main -p $NUM_PROCS -superlu < "${RB_DIR}"/"${FILE}" > /dev/null
        echo >> $OUT_FILE
        # break
    done
done

cd ../..