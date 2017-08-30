all: benchmark_count benchmark_locate benchmark_extract

benchmark_count: benchmark_count.cpp
	g++ -std=c++11 -O3 -DNDEBUG -I ~/include -L ~/lib benchmark_count.cpp -o benchmark_count -lsdsl -ldivsufsort -ldivsufsort64

benchmark_locate: benchmark_locate.cpp
	g++ -std=c++11 -O3 -DNDEBUG -I ~/include -L ~/lib benchmark_locate.cpp -o benchmark_locate -lsdsl -ldivsufsort -ldivsufsort64

benchmark_extract: benchmark_extract.cpp
	g++ -std=c++11 -O3 -DNDEBUG -I ~/include -L ~/lib benchmark_extract.cpp -o benchmark_extract -lsdsl -ldivsufsort -ldivsufsort64
