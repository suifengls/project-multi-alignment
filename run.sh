#!/bin/bash

echo "----- Test case one -----"
./multi test1/A.txt test1/B.txt test1/C.txt  
# max 5262 
./multi test2/NM_000558.txt test2/NM_008218.txt test2/NM_013096.txt
#echo "----- Test case two -----"
# max 12362
./multi test3/NM_001030004.txt test3/NM_178850.2.txt test3/XM_514664.txt
#echo "----- Test case three -----"
#max 8793
#./threeseq test4/NM_001243563.1.txt test4/NM_010019.3.txt test4/NM_014326.3.txt 

#echo "----- Test case four -----"
#./threeseq test5/NM_000457.4.txt test5/NM_000545.5.txt test5/NM_008261.2.txt

#echo "----- Test case five -----"
#./threeseq test6/NM_000492.3.txt test6/NM_021050.2.txt test6/NM_031506.1.txt
