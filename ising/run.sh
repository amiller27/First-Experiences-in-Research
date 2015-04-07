g++ ising_2d.cc -o ising -pthread --std=c++11 && \
./ising && \
python -c "import ising_plotter; ising_plotter.plot(\"ising_2d.dat\");"

