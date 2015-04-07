g++ ising_3d.cc -o ising -pthread --std=c++11 && \
./ising && \
python -c "import ising_plotter; ising_plotter.plot(\"ising_3d.dat\");"
