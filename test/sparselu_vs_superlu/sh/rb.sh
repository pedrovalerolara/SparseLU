TEST_DIR=$HOME/workplan-summer-2021-internship-copy/LASs-sparse-npgetrf-and-dense-tile-pivoting/LASs-sparse-npgetrf/test/slass_v_superlu
RB_DIR=$TEST_DIR/matrices/Rutherford_Boeing

cd $RB_DIR

for file in ./*.rb; do
    $TEST_DIR/matrix_attributes $file
done

cd $TEST_DIR