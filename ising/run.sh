g++ ising_2d.cc -o ising -pthread --std=c++11 && \
./ising && \
python -c "import ising_2d; ising_2d.plot(\"ising_2d.dat\");"

