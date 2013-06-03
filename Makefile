all:
		g++ -O3 -I ${BOOST_ROOT} multi.cpp -o multi
clean:
		rm -rf multi
