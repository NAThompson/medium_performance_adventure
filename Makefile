all: bench.x

bench.x: bench.cpp
	clang++ -g -fno-omit-frame-pointer -fno-unroll-loops --std=c++17 -march=native -O3 -Wfatal-errors -I../boost  -I/usr/local/include bench.cpp -L/usr/local/lib -lbenchmark -lbenchmark_main -pthread
