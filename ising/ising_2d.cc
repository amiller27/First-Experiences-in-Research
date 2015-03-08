#include <cmath>
#include <iostream>
#include <sstream>
#include <random>
#include <fstream>
#include <thread>
#include <future>

// std::ifstream randfile ("random.txt");

// set up random number aliases
std::default_random_engine generator;
std::uniform_real_distribution<double> distribution(0.0, 1.0);

// holds data passed to and returned from swapTemps function
struct WalkData {
	double T;
	int U_total;
	int magnetization_total;
};

// Simulator class declaration
class Simulator
{
	// BASIC SIMULATION DATA
	int steps_per_walk;
	int steps_between_swaps;
	int max_flips;
	int lattice_size;
	double* k_T_values;
	int num_of_simulations;

	// each item is either -1 or 1
	// 1 if next swap should be with the one above, -1 if next swap should be with the one below
	// 
	// Note: the temperatures on the ends do not swap if their swap value in this array points toward the outside
	int* nextSwaps;

	// arrays of the current data at each temperature, updated by swapTemps
	int* energies;
	WalkData* current_walk_data;

	// array of num_of_simulations-1 random numbers, generated for each swap
	// 
	// these must be stored so that the two sides of the swap
	// for the same two temperatures use the same random number
	// and get the same result
	double* randoms;

	// FLAGS USED FOR SYNCHRONIZATION WHILE SWAPPING TEMPERATURES
	
	// the number of threads which are ready to swap, or the number of threads which are finished swapping
	std::atomic<int> ready;

	// variables used with a condition_variable for signalling
	std::mutex mtx;
	std::condition_variable cv;
	bool done;
	bool go;

	// flag indicating whether a master thread has already been picked for the swap
	std::atomic_flag master_waiting ATOMIC_FLAG_INIT;

public:
	Simulator(int steps_per_walk, int steps_between_swaps, int max_flips, int lattice_size, double k_T_values[], int num_of_simulations);
	double** simulate();
	WalkData swapTemps(int U, WalkData data);
};

void printLattice(const int* const* const* lattice, int lattice_size) {
	for (int i = 1; i < lattice_size + 1; i++) {
		for (int j = 1; j < lattice_size + 1; j++) {
			std::cout << (*lattice[i][j] == -1 ? " -1 " : std::string("  ") + std::to_string(*lattice[i][j]) + " ");
		}
		std::cout << std::endl;
	}
}

void walk(int steps, int steps_between_swaps, int flips_per_iteration, double k_times_T, int lattice_size, std::promise<double*>* result_promise, Simulator* sim) {

	// 3D array of int*s, initially all up
	// 1 is up, -1 is down
	// note: all edges are pointers to the same ints represented by the opposite edge
	// corners are NOT initialized
	int*** lattice = new int**[lattice_size + 2];
	for (int i = 1; i < lattice_size + 1; i++) {
		lattice[i] = new int*[lattice_size + 2];
		for (int j = 1; j < lattice_size + 1; j++) {
			lattice[i][j] = new int(1);
		}

		// initialize edge cases
		lattice[i][0] = lattice[i][lattice_size];
		lattice[i][lattice_size + 1] = lattice[i][1];
	}

	// initialize edge cases
	lattice[0] = new int*[lattice_size + 2];
	lattice[lattice_size + 1] = new int*[lattice_size + 2];
	for (int j = 1; j < lattice_size + 1; j++) {
		lattice[0][j] = lattice[lattice_size][j];
		lattice[lattice_size + 1][j] = lattice[1][j];
	}

	// printLattice(lattice, lattice_size);

	// potential
	int U = -4 * lattice_size * lattice_size;
	int U_total = 0;
	// number of up-spin particles
	int num_up = lattice_size * lattice_size;
	// magnetization
	int magnetization = lattice_size * lattice_size;
	int magnetization_total = 0;

	// counter to measure when to flip temperatures
	int current_run_counter = 0;

	// next lattice point to flip
	int current_point = 0;
	for (int step = 0; step < steps; step++, current_run_counter++) {

		// std::stringstream ss;
		// ss << "[Thread " << std::this_thread::get_id() << "]\twalking " << step << std::endl;
		// std::cout << ss.str();

		if (current_run_counter == steps_between_swaps) {
			WalkData old_data;
			old_data.T = k_times_T;
			old_data.U_total = U_total;
			old_data.magnetization_total = magnetization_total;

			WalkData new_data = sim->swapTemps(U, old_data);

			k_times_T = new_data.T;
			U_total = new_data.U_total;
			magnetization_total = new_data.magnetization_total;

			current_run_counter = 0;
		}

		// calculate change in values from flips
		int delta_U = 0;
		int delta_num_up = 0;
		for (int flip = 0; flip < flips_per_iteration; flip++, current_point = (current_point + 1)%(lattice_size * lattice_size)) {
			int current_x = current_point / lattice_size;
			int current_y = current_point % lattice_size;
			// if they're in the same state initially, the energy goes up
			// otherwise, it goes down
			delta_U += (*lattice[current_x + 1][current_y + 1] == *lattice[current_x + 1][current_y + 1 + 1] ? 2 : -2);
			delta_U += (*lattice[current_x + 1][current_y + 1] == *lattice[current_x + 1][current_y - 1 + 1] ? 2 : -2);
			delta_U += (*lattice[current_x + 1][current_y + 1] == *lattice[current_x + 1 + 1][current_y + 1] ? 2 : -2);
			delta_U += (*lattice[current_x + 1][current_y + 1] == *lattice[current_x - 1 + 1][current_y + 1] ? 2 : -2);

			// debugging
			// std::cout << current_point << std::endl;
			// std::cout << current_x << std::endl;
			// std::cout << current_y << std::endl;
			// flip the particle
			*lattice[current_x + 1][current_y + 1] = -1 * (*lattice[current_x + 1][current_y + 1]);
			delta_num_up += *lattice[current_x + 1][current_y + 1];
		}

		// decide whether to accept
		// std::cout << delta_U << " " << k_times_T << std::endl;
		// std::string line;
		// getline(randfile, line);
		// double nextrand = atof(line.substr(0, line.length()).c_str());
		double nextrand = distribution(generator);
		// std::cout << nextrand << " " << exp(-delta_U / k_times_T) << std::endl;
		if (delta_U <= 0 || nextrand < exp(-delta_U / k_times_T)) {
			// accept
			// std::cout << "accepted" << std::endl;
			U += delta_U;
			num_up += delta_num_up;
			magnetization = num_up - (lattice_size * lattice_size - num_up);
		} else {
			// reject
			// std::cout << "rejected" << std::endl;
			for (int flip = (lattice_size * lattice_size + current_point - flips_per_iteration)%(lattice_size * lattice_size); flip != current_point; flip = (flip + 1)%(lattice_size * lattice_size)) {
				int current_x = flip / lattice_size;
				int current_y = flip % lattice_size;
				// debugging
				// std::cout << flip << std::endl;
				// std::cout << current_x << std::endl;
				// std::cout << current_y << std::endl;
				*lattice[current_x + 1][current_y + 1] = -1 * (*lattice[current_x + 1][current_y + 1]);
			}
		}

		// debugging
		// printLattice(lattice, lattice_size);
		// std::cout << num_up << " " << U << " " << magnetization << std::endl;

		// add values into running totals
		U_total += U;
		magnetization_total += abs(magnetization);

		// print for debugging
		// printLattice(lattice, lattice_size);
	}

	// printLattice(lattice, lattice_size);

	//clean up
	delete[] lattice[0];
	for (int i = 1; i < lattice_size + 1; i++) {
		for (int j = 1; j < lattice_size + 1; j++) {
			delete lattice[i][j];
		}
		delete[] lattice[i];
	}
	delete[] lattice[lattice_size + 1];
	delete[] lattice;

	double* result_array = new double[3] {k_times_T, double(U_total) / steps, double(magnetization_total) / steps};
	(*result_promise).set_value(result_array);
}

Simulator::Simulator (
	int steps_per_walk,
	int steps_between_swaps,
	int max_flips,
	int lattice_size,
	double k_T_values[],
	int num_of_simulations
) : steps_per_walk(steps_per_walk),
	steps_between_swaps(steps_between_swaps),
	max_flips(max_flips),
	lattice_size(lattice_size),
	k_T_values(k_T_values),
	num_of_simulations(num_of_simulations)
{
	nextSwaps = new int[num_of_simulations];
	for (int i = 0; i < num_of_simulations; i += 2) {
		nextSwaps[i] = 1;
	}
	for (int i = 1; i < num_of_simulations; i += 2) {
		nextSwaps[i] = -1;
	}

	ready = 0;

	done = false;
	go = false;

	randoms = new double[num_of_simulations / 2];

	energies = new int[num_of_simulations];
	for (int i = 0; i < num_of_simulations; i++) {
		energies[i] = 0;
	}

	current_walk_data = new WalkData[num_of_simulations];
}

double** Simulator::simulate() {

	double** results = new double*[3];
	results[0] = new double[num_of_simulations];
	results[1] = new double[num_of_simulations];
	results[2] = new double[num_of_simulations];

	std::thread** active_threads = new std::thread*[num_of_simulations];
	std::promise<double*>** active_promises = new std::promise<double*>*[num_of_simulations];

	// start threads
	for (int i = 0; i < num_of_simulations; i++) {
		active_promises[i] = new std::promise<double*>();
		active_threads[i] = new std::thread(walk, steps_per_walk, steps_between_swaps, max_flips, k_T_values[i], lattice_size, active_promises[i], this);
	}

	for (int i = 0; i < num_of_simulations; i++) {

		// store results
		double* next_result = (*active_promises[i]).get_future().get();

		int next_index = -1;
		// maybe change this to binary search in the future
		// for now, temperature sets are small enough to use sequential
		for (int i = 0; i < num_of_simulations; i++) {
			if (k_T_values[i] == next_result[0]) {
				next_index = i;
				break;
			}
		}

		results[0][next_index] = next_result[0];
		results[1][next_index] = next_result[1];
		results[2][next_index] = next_result[2];
		// std::cout << next_index << std::endl;
		// std::cout << pending_futures << std::endl;

		// plot them
		std::stringstream ss;
		ss << "Temperature parameter: " << next_result[0] << std::endl;
		ss << "Average U:\t\t\t" << next_result[1] << std::endl;
		ss << "Average Magnetization:\t\t\t\t" << next_result[2] << std::endl;
		std::cout << ss.str();
		delete[] next_result;
	}

	for (int i = 0; i < num_of_simulations; i++) {
		active_threads[i]->join();
		delete active_threads[i];
		delete active_promises[i];
	}
	delete[] active_threads;
	delete[] active_promises;

	return results;
}

WalkData Simulator::swapTemps(int U, WalkData data) {
	double T = data.T;

	int temperatureIndex = -1;
	// maybe change this to binary search in the future
	// for now, temperature sets are small enough to use sequential
	for (int i = 0; i < num_of_simulations; i++) {
		if (k_T_values[i] == T) {
			temperatureIndex = i;
			break;
		}
	}

	// std::stringstream ss;
	// ss << "[Thread " << std::this_thread::get_id() << "]\tWaiting to swap " << U << " " << T << " " << nextSwaps[temperatureIndex] << std::endl;
	// std::cout << ss.str();

	int otherTemperatureIndex = temperatureIndex + nextSwaps[temperatureIndex];
	nextSwaps[temperatureIndex] *= -1;

	energies[temperatureIndex] = U;
	current_walk_data[temperatureIndex] = data;

	bool isMaster = false;

	ready++;
	if (!master_waiting.test_and_set()) {
		// this is the first thread to finish
		// make this thread the master thread
		// std::stringstream ss1;
		// ss1 << "[Thread " << std::this_thread::get_id() << "]\tMaster thread " << U << " " << T << " " << nextSwaps[temperatureIndex] << std::endl;
		// std::cout << ss1.str();
		
		isMaster = true;

		while (ready < num_of_simulations) {
			// std::stringstream ss2;
			// ss2 << "[Thread " << std::this_thread::get_id() << "]\tWaiting for ready " << ready  << " " << U << " " << T << " " << nextSwaps[temperatureIndex] << std::endl;
			// std::cout << ss2.str();
		}

		for (int i = 0; i < num_of_simulations / 2; i++) {
			randoms[i] = distribution(generator);
		}

		std::unique_lock<std::mutex> lck (mtx);
		done = false;
		ready = 0;
		go = true;
		cv.notify_all();

	} else {
		// std::stringstream ss1;
		// ss1 << "[Thread " << std::this_thread::get_id() << "]\tNormal thread " << ready  << " " << U << " " << T << " " << nextSwaps[temperatureIndex] << std::endl;
		// std::cout << ss1.str();
		std::unique_lock<std::mutex> lck (mtx);
		while (!go) {
			// std::stringstream ss1;
			// ss1 << "[Thread " << std::this_thread::get_id() << "]\tNormal thread waiting for ready " << ready  << " " << U << " " << T << " " << nextSwaps[temperatureIndex] << std::endl;
			// std::cout << ss1.str();
			cv.wait(lck);
		}
	}

	WalkData result = data;
	if (otherTemperatureIndex == -1 || otherTemperatureIndex == num_of_simulations) {
		// std::stringstream ss14;
		// ss14 << "[Thread " << std::this_thread::get_id() << "]\tReturning without swapping" << std::endl;
		// std::cout << ss14.str();
	} else {

		int otherEnergy = energies[otherTemperatureIndex];
		int dU = U - otherEnergy;
		double dTReciprocal = 1/T - 1/k_T_values[otherTemperatureIndex];
		double probabilityOfSwap = fmin(1.0, exp(dU * dTReciprocal));
		double random = randoms[(temperatureIndex < otherTemperatureIndex ? temperatureIndex : otherTemperatureIndex) / 2];
		if (random < probabilityOfSwap) {
			result = current_walk_data[otherTemperatureIndex];
		}
	}

	ready++;

	if (isMaster) {
		while (ready < num_of_simulations) {}
		ready = 0;

		std::unique_lock<std::mutex> lck (mtx);
		master_waiting.clear();
		done = true;
		go = false;
		cv.notify_all();
	} else {
		std::unique_lock<std::mutex> lck (mtx);
		while (!done) {
			cv.wait(lck);
		}
	}

	// std::stringstream ss8;
	// ss8 << "[Thread " << std::this_thread::get_id() << "]\tDone swapping " << U << " " << T << " " << result << std::endl;
	// std::cout << ss8.str();

	return result;
}

int main(int argc, char const *argv[]) {

	// constants
	int num_of_temps = 100;
	double low = 0.1;
	double increment = 0.1;

	int steps_per_walk = 1000000;
	int steps_between_swaps = 100;
	int max_flips = 1;
	int lattice_size = 4;
	
	// set up temperature arguments
	double temps[num_of_temps];
	for (int i = 0; i < num_of_temps; i++) {
		temps[i] = low + i*increment;
	}

	// simulate
	Simulator* simulator = new Simulator(steps_per_walk, steps_between_swaps, max_flips, lattice_size, temps, num_of_temps);
	double** results = simulator->simulate();
	delete simulator;

	// write data to file
	std::ofstream data_file ("ising_2d.dat");
	for (int i = 0; i < num_of_temps; i++) {
		data_file << temps[i] << " " << results[1][i] << " " << results[2][i] << std::endl;
	}
	data_file.close();

	// clean up
	delete[] results[0];
	delete[] results[1];
	delete[] results;

	return 0;
}