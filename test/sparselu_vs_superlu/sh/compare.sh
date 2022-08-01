#dgterf

# TEST ALL MATRICES AND OUTPUT IN TAB-DELIMITED
NUM_PROCS=30
# SIZES=(1024 2048 4096 8192 16384 20480 25600 32768)
# SIZES=(1024 2048 3072 4096 5120 6144 7168 8192 9216 10240 11264 12288 13312 14336 15360 16384 17480 18432 19456 20480 21504 22528 23552 24576 25600 26624 27648 28672 29696 30720 31744 32768)
SIZES=(17408)
DENSITIES=(0.001 0.0015 0.002 0.0025 0.003 0.0035 0.004 0.0045 0.005 0.0055 0.006 0.0065 0.007 0.0075)
LIBRARIES=(lass superlu)
TEST_DIR=$HOME/workplan-summer-2021-internship-copy/LASs-sparse-npgetrf-and-dense-tile-pivoting/LASs-sparse-npgetrf/test/slass_v_superlu
OUT_FILE=$TEST_DIR/out/compare-lass-only.txt
RB_DIR=$TEST_DIR/matrices/Synthetic

# rm $OUT_FILE
# echo -e 'File\tSize\tDensity\tLibrary\tRuntime' >> $OUT_FILE

cd "${RB_DIR}"
for LIBRARY in "${LIBRARIES[@]}"; do
    for DENSITY in "${DENSITIES[@]}"; do
        for SIZE in "${SIZES[@]}"; do
            FILE="${SIZE}-${DENSITY}".rb

            # Size
            echo -e -n "${SIZE}\t">> $OUT_FILE

            # Density
            echo -e -n "$DENSITY\t">> $OUT_FILE

            # Library
            echo -e -n "${LIBRARY}\t">> $OUT_FILE

            # Runtime
            numactl --interleave=all $TEST_DIR/main -p $NUM_PROCS -${LIBRARY} < "${RB_DIR}"/"${FILE}" > /dev/null &
            wait
            echo >> $OUT_FILE
        done
    done
done
# for LIBRARY in "${LIBRARIES[@]}"; do
#     for FILE in ./*.rb; do
#         # Matrix attributes
#         echo $FILE
        
#         $TEST_DIR/matrix_attributes $FILE
#         echo -e -n "\t" >> $OUT_FILE

#         # Library
#         echo -e -n "${LIBRARY}\t">> $OUT_FILE

#         # Runtime
#         numactl --interleave=all $TEST_DIR/main -p $NUM_PROCS -${LIBRARY} < "${RB_DIR}"/"${FILE}" > /dev/null
#         echo >> $OUT_FILE
#     done
# done

cd ../..
