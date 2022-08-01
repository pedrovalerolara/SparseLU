rm ./out.txt

echo 1000x1000 >> out.txt
./test_sparse_dnpgetrf 1000 1000 &>> out.txt 
echo "finished 1000x1000"
echo "/************************************************/" >> out.txt
echo 2000x2000 >> out.txt
./test_sparse_dnpgetrf 2000 2000 &>> out.txt
echo "finished 2000x2000"
echo "/************************************************/" >> out.txt
echo 4000x4000 >> out.txt
./test_sparse_dnpgetrf 4000 4000 &>> out.txt 
echo "finished 4000x4000"
echo "/************************************************/" >> out.txt
echo 8000x8000 >> out.txt
./test_sparse_dnpgetrf 8000 8000 &>> out.txt 
echo "finished 8000x8000"
echo "/************************************************/" >> out.txt
echo 16000x16000 >> out.txt
./test_sparse_dnpgetrf 16000 16000 &>> out.txt 
echo "finished 16000x16000"
echo "/************************************************/" >> out.txt
echo 32000x32000 >> out.txt
./test_sparse_dnpgetrf 32000 32000 &>> out.txt
echo "finished 32000x32000"
echo "/************************************************/" >> out.txt
echo 64000x64000 >> out.txt
./test_sparse_dnpgetrf 64000 64000 &>> out.txt 
echo "finished 64000x64000"
