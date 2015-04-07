import numpy as np
import math
import Tkinter
import threading

# random_file = open("random.txt")

class App(object):
	def __init__(self, master, lattice_size):
		self.size = 512

		self.main_canvas = Tkinter.Canvas(master, width = self.size, height = self.size)
		self.main_canvas.pack()

		self.side_length = self.size / lattice_size
		self.lattice = []
		for i in range(lattice_size):
			self.lattice.append([])
			for j in range(lattice_size):
				x = self.side_length * i
				y = self.side_length * j
				self.lattice[i].append(self.main_canvas.create_rectangle(x, y, x + self.side_length, y + self.side_length, fill = "blue"))

	def draw(self, square_location, spin):
		self.main_canvas.itemconfig(self.lattice[square_location[0]][square_location[1]], fill = ("blue" if spin == 1 else "red"))

class SimulationThread(threading.Thread):
	def __init__(self, result, k_T_values, lattice_size, gui = None):
		super(SimulationThread, self).__init__()

		self.result = result

		self.k_T_values = k_T_values
		self.lattice_size = lattice_size
		self.gui = gui

	def run(self):
		self.result[0], self.result[1] = simulate(self.k_T_values, self.lattice_size, self.gui)

def printLattice(lattice):
	for row in lattice:
		for item in row:
			print (" -1" if item==-1 else "  1"),
		print

def walk(steps, max_flips_per_iteration, k_times_T, lattice_size, plotter = None):
	lattice = list(list(1 for j in range(lattice_size)) for i in range(lattice_size))

	U = -4 * lattice_size**2
	U_total = 0

	num_up = lattice_size**2

	magnetization = lattice_size**2
	magnetization_total = 0

	current_point = 0
	for step in xrange(steps):
		delta_U = 0
		delta_num_up = 0
		flips_per_iteration = np.random.random_integers(1, max_flips_per_iteration)
		for flip in range(flips_per_iteration):
			current_x = current_point / lattice_size
			current_y = current_point % lattice_size

			if current_x != lattice_size - 1:
				delta_U += 2 if lattice[current_x][current_y]==lattice[current_x + 1][current_y] else -2
			else:
				delta_U += 2 if lattice[current_x][current_y]==lattice[0][current_y] else -2

			delta_U += 2 if lattice[current_x][current_y]==lattice[current_x - 1][current_y] else -2

			if current_y != lattice_size - 1:
				delta_U += 2 if lattice[current_x][current_y]==lattice[current_x][current_y + 1] else -2
			else:
				delta_U += 2 if lattice[current_x][current_y]==lattice[current_x][0] else -2

			delta_U += 2 if lattice[current_x][current_y]==lattice[current_x][current_y - 1] else -2

			lattice[current_x][current_y] = -lattice[current_x][current_y]
			delta_num_up += lattice[current_x][current_y]

			current_point = (current_point + 1) % lattice_size**2

		# rand = float(random_file.next()[:-2])
		rand = np.random.random()
		# print delta_U, k_times_T
		# print rand, math.exp(-float(delta_U) / k_times_T)
		if delta_U <= 0 or rand < math.exp(-float(delta_U) / k_times_T):
			U += delta_U
			num_up += delta_num_up
			magnetization = 2 * num_up - lattice_size**2
		else:
			flip = (current_point - flips_per_iteration) % (lattice_size**2)
			while flip != current_point:
				current_x = flip / lattice_size
				current_y = flip % lattice_size
				lattice[current_x][current_y] = -lattice[current_x][current_y]
				flip = (flip + 1) % (lattice_size**2)
		
		U_total += U
		magnetization_total += abs(magnetization)

		# for debugging
		# printLattice(lattice)
		# print num_up, U, magnetization
		# print
		if plotter:
			last_x = (current_point - 1)%(lattice_size**2) / lattice_size
			last_y = (current_point - 1)%(lattice_size**2) % lattice_size
			plotter.draw([last_x, last_y], lattice[last_x][last_y])
			# print step

	return (float(U_total) / steps, float(magnetization_total) / steps)

def simulate(k_T_values, lattice_size, gui = None):
	energies = []
	magnetizations = []
	for k_T in k_T_values:
		energy, magnetization = walk(100000, 5, k_T, lattice_size, gui)
		print "Temperature parameter: %f"%k_T
		print "Average U:\t\t\t%f"%energy
		print "Average Magnetization:\t\t\t\t%f"%magnetization
		energies.append(energy)
		magnetizations.append(magnetization)
	return energies, magnetizations

def temp_range(low, high, increment = 1):
	temps = []
	while low < high:
		temps.append(low)
		low += increment
	return temps

if __name__=="__main__":
	lattice_size = 16
	temps = temp_range(0.2, 4, 0.2)

	root = Tkinter.Tk()
	app = App(root, lattice_size)

	distributions = [[], []]

	sim_thread = SimulationThread(distributions, temps, lattice_size, app)
	sim_thread.start()

	root.mainloop()

	energy_distribution, magnetization_distribution = distributions

	data_file = open("ising_2d.dat", 'w')
	for i in range(len(temps)):
		data_file.write("%f %f %f\n"%(temps[i], energy_distribution[i], magnetization_distribution[i]))
	data_file.close()

	import ising_plotter
	ising_plotter.plot("ising_2d.dat")
