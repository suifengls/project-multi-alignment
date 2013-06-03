all:
		g++ -O2 -I ${BOOST_ROOT} threeseq.cpp -o threeseq
clean:
		rm -rf threeseq
