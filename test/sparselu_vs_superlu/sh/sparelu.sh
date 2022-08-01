#dgterf

# TEST ALL MATRICES AND OUTPUT IN TAB-DELIMITED
SIZES=(1024 2048 3072 4096 5120 6144 7168 8192 9216 10240 11264 12288 13312 14336 15360)
DENSITY=0.005
TILE_SIZE=1024
TEST_DIR=$HOME/workplan-summer-2021-internship-copy/LASs-sparse-npgetrf-and-dense-tile-pivoting/LASs-sparse-npgetrf/test/slass_v_superlu
OUT_FILE=$TEST_DIR/out/lass_results.txt
RB_DIR=$TEST_DIR/matrices/Synthetic
INDEX=60

# rm $OUT_FILE
# echo -e 'INDEX\tFILE\tSIZE\tDENSITY\tTILE_SIZE\tRUNTIME_SLASS' >> $OUT_FILE

cd "${RB_DIR}"
for SIZE in "${SIZES[@]}"; do
    FILE="${SIZE}".rb

    ((INDEX=INDEX+1))
    
    # Index
    echo -e -n $INDEX'\t'>> $OUT_FILE

    # File
    echo -e -n $FILE'\t'>> $OUT_FILE

    # Matrix Attributes
    $TEST_DIR/matrix_attributes "${RB_DIR}"/"${FILE}"
    echo -e -n '\t'>> $OUT_FILE

    # Size
    echo -e -n "${SIZE}\t" >> $OUT_FILE

    # Density
    echo -e -n $DENSITY>> $OUT_FILE
    echo -e -n '\t'>> $OUT_FILE

    # Tile Size
    echo -e -n $TILE_SIZE>> $OUT_FILE
    echo -e -n '\t'>> $OUT_FILE

    # Test Run
    $TEST_DIR/main -p $NUM_PROCS -slass < "${RB_DIR}"/"${FILE}" > /dev/null
    echo >> $OUT_FILE
    # break
done

cd ../..