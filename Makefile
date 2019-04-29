all: bench.x

bench.x: bench.cpp
	 g++-8 --std=c++17 -march=native -O3 -Wfatal-errors -I../boost  -I/usr/local/include bench.cpp -L/usr/local/lib -lbenchmark -lbenchmark_main -pthread
