# Generate variously-sized matrices 
# SIZES=(1024 2048 4096 8192 16384 20480 25600 32768)
#SIZES=(1024 2048 3072 4096 5120 6144 7168 8192 9216 10240 11264 12288 13312 14336 15360 16384 17480 18432 19456 20480 21504 22528 23552 24576 25600 26624 27648 28672 29696 30720 31744 32768)
SIZES=(100)
# DENSITIES=(0.001 0.0015 0.002 0.0025 0.003 0.0035 0.004 0.0045 0.005 0.0055 0.006 0.0065 0.007 0.0075)
DENSITIES=(0.05)
TEST_DIR=$HOME/workplan-summer-2021-internship-copy/LASs-sparse-npgetrf-and-dense-tile-pivoting/LASs-sparse-npgetrf/test/slass_v_superlu
RB_DIR=$TEST_DIR/matrices/Synthetic

cd $RB_DIR

# Generate the synthetic matrices
for SIZE in "${SIZES[@]}"; do
	for DENSITY in "${DENSITIES[@]}"; do
		echo "${SIZE}-${DENSITY}"
		$TEST_DIR/matrix_gen -mn "${SIZE}" -d "${DENSITY}" -of "${RB_DIR}"/"${SIZE}-${DENSITY}".rb
	done
done

cd ../..
