import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.pylab as plab
import threading

import ui
import scene

# random_file = open("random.txt")

class LatticeScene(scene.Scene):
	def __init__(self, lattice_size):
		super(LatticeScene, self).__init__()
		self.lattice_size = lattice_size
		self.updates = []

	def setup(self):
		self.spacing = (min(self.size.w, self.size.h)) / self.lattice_size
		scene.stroke(1, 1, 1)
		scene.stroke_weight(1)
		scene.fill(0, 0, 1)
		for i in range(self.lattice_size):
			for j in range(self.lattice_size):
				scene.rect(i*self.spacing, j*self.spacing, self.spacing, self.spacing)

	def update_square(self, location, spin):
		# print "updating", location, spin
		self.updates.append((location, spin))

	def reset(self, lattice_size):
		self.lattice_size = lattice_size
		self.updates = []
		self.setup()

	def draw(self):
		scene.stroke(1, 1, 1)
		scene.stroke_weight(1)

		for location, spin in self.updates:
			if spin == 1:
				scene.fill(0, 0, 1)
			else:
				scene.fill(1, 0, 0)
			scene.rect(location[0]*self.spacing, location[1]*self.spacing, self.spacing, self.spacing)
		self.updates = []

class PlotterScene(scene.Scene):
	def __init__(self,):
		super(PlotterScene, self).__init__()
		self.first_frame_drawn = False
		self.first_point_drawn = False
		self.point_queue = []

	def reset(self, min_temp, max_temp, max_magnetization, max_energy):
		self.min_temp = min_temp
		self.max_temp = max_temp
		self.max_magnetization = max_magnetization
		self.max_energy = max_energy

		self.first_frame_drawn = False
		self.first_point_drawn = False
		self.point_queue = []
		self.last_point = ()

	def add_point(self, temperature, magnetization, energy):
		print 'adding'
		self.point_queue.append((temperature, magnetization, energy))
		print self.point_queue

	def draw(self):
		if not self.first_frame_drawn:
			scene.background(1, 1, 1)

			scene.stroke(0, 0, 0)
			scene.stroke_weight(2)

			# lower graph axes
			scene.line(0, self.size.h / 4 - 4, self.size.w, self.size.h / 4 - 4)
			scene.line(1, 0, 1, self.size.h / 2 - 8)

			# upper graph axes
			scene.line(0, self.size.h / 2 + 8, self.size.w, self.size.h / 2 + 8)
			scene.line(1, self.size.h / 2 + 8, 1, self.size.h)

			self.first_frame_drawn = True

		scene.stroke(0, 1, 1)
		scene.stroke_weight(4)
		scene.fill(0, 1, 1)

		if not self.first_point_drawn and self.point_queue:
			self.last_point = self.point_queue.pop(0)
			x = (self.last_point[0] - self.min_temp) * self.size.w / (self.max_temp - self.min_temp)
			y_magnetization = (self.size.h / 2 - 8) * self.last_point[1] / self.max_magnetization + (self.size.h / 2 + 8)
			y_energy = (self.size.h / 4 - 4) * self.last_point[2] / self.max_energy + (self.size.h / 4 - 4)
			scene.ellipse(x, y_magnetization, 2, 2)
			scene.ellipse(x, y_energy, 2, 2)

			self.first_point_drawn = True

		for temperature, magnetization, energy in self.point_queue:
			x1 = (self.last_point[0] - self.min_temp) * self.size.w / (self.max_temp - self.min_temp)
			x2 = (temperature - self.min_temp) * self.size.w / (self.max_temp - self.min_temp)
			y1_magnetization = (self.size.h / 2 - 8) * self.last_point[1] / self.max_magnetization + (self.size.h / 2 + 8)
			y2_magnetization = (self.size.h / 2 - 8) * magnetization / self.max_magnetization + (self.size.h / 2 + 8)
			y1_energy = (self.size.h / 4 - 4) * self.last_point[2] / self.max_energy + (self.size.h / 4 - 4)
			y2_energy = (self.size.h / 4 - 4) * energy / self.max_energy + (self.size.h / 4 - 4)
			scene.line(x1, y1_magnetization, x2, y2_magnetization)
			scene.line(x1, y1_energy, x2, y2_energy)

			print x1, x2, y1_magnetization, y2_magnetization, y1_energy, y2_energy

			self.last_point = (temperature, magnetization, energy)
		self.point_queue = []

class SimulationThread(threading.Thread):
	def __init__(self, result, k_T_values, lattice_size, gui = None, plotter = None, on_temperature_changed = None):
		super(SimulationThread, self).__init__()

		self.result = result

		self.k_T_values = k_T_values
		self.lattice_size = lattice_size
		self.gui = gui
		self.plotter = plotter
		self.on_temperature_changed = on_temperature_changed

		self.should_stop = threading.Event()

	def run(self):
		energy_distribution, magnetization_distribution = simulate(self.k_T_values, self.lattice_size, self.gui, self.plotter, on_temperature_changed)

		data_file = open("ising_2d.dat", 'w')
		for i in range(len(self.k_T_values)):
			data_file.write("%f %f %f\n"%(self.k_T_values[i], energy_distribution[i], magnetization_distribution[i]))
		data_file.close()

	def stop(self):
		self.should_stop.set()

def printLattice(lattice):
	for row in lattice:
		for item in row:
			print (" -1" if item==-1 else "  1"),
		print

def walk(steps, max_flips_per_iteration, k_times_T, lattice_size, gui = None):
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
		if gui:
			last_x = (current_point - 1)%(lattice_size**2) / lattice_size
			last_y = (current_point - 1)%(lattice_size**2) % lattice_size
			gui.update_square([last_x, last_y], lattice[last_x][last_y])
			# print step

	return (float(U_total) / steps, float(magnetization_total) / steps)

def simulate(k_T_values, lattice_size, gui = None, plotter = None, on_temperature_changed = None):
	energies = []
	magnetizations = []
	for k_T in k_T_values:
		try:
			on_temperature_changed(k_T)
		except TypeError:
			# the function was not passed in, so do nothing
			pass
		energy, magnetization = walk(30000, 1, k_T, lattice_size, gui)
		print "Temperature parameter: %f"%k_T
		print "Average U:\t\t\t%f"%energy
		print "Average Magnetization:\t\t\t\t%f"%magnetization
		energies.append(energy)
		magnetizations.append(magnetization)

		if plotter:
			plotter.add_point(k_T, magnetization, energy)
	return energies, magnetizations

def plot(filename):
	data_file = open(filename)
	data = map(str.split, data_file.read().split("\n"), " ")[:-1]
	data_file.close()
	temps = list(float(i[0]) for i in data)
	energies = list(float(i[1]) for i in data)
	magnetizations = list(float(i[2]) for i in data)


	plt.subplot(2, 1, 2)
	plt.plot(temps, energies)
	plt.ylabel("Energy")
	plt.subplot(2, 1, 1)
	plt.plot(temps, magnetizations)
	plt.ylabel("Magnetization")
	plt.xlabel("Temperature (k * T / epsilon)")
	plt.show()

def temp_range(low, high, increment = 1):
	temps = []
	while low <= high:
		temps.append(low)
		low += increment
	return temps

if __name__=="__main__":
	def on_temperature_changed(temperature):
		main_view['running_view']['temperature_label'].text = str(temperature)

	def on_lattice_size_changed(sender):
		main_view['lattice_size_label'].text = str(2**int(sender.value*4 + 1))
		lattice_scene.reset(2**int(sender.value*4 + 1))

	def on_reset_pressed(sender):
		sim_thread.stop()
		sim_thread.join()
		main_view['running_view'].send_to_back()
		lattice_scene.reset(2**int(main_view['lattice_size_slider'].value * 4 - 1))

	def on_start_pressed(sender):
		try:
			low_temp = float(main_view['initial_temperature_field'].text)
			high_temp = float(main_view['final_temperature_field'].text)
		except ValueError:
			main_view['temperature_error'].bring_to_front()
			return
		if high_temp == low_temp:
			main_view['temperature_error'].bring_to_front()
			return
		if high_temp < low_temp:
			t = high_temp
			high_temp = low_temp
			low_temp = t

		lattice_size = 2**int(main_view['lattice_size_slider'].value * 4 + 1)

		main_view['temperature_error'].send_to_back()
		temps = temp_range(low_temp, high_temp, (high_temp - low_temp) / 12)

		lattice_scene.reset(lattice_size)
		plotter_scene.reset(temps[0], temps[-1], lattice_size**2, 4 * lattice_size**2)

		main_view['running_view'].bring_to_front()

		distributions = [[], []]
		# NEED TO FIGURE THIS OUT WITHOUT GLOBAL VARIABLE
		# MAYBE THROW EVERYTHING IN A CLASS
		global sim_thread
		sim_thread = SimulationThread(distributions, temps, lattice_size, lattice_scene, plotter_scene)
		sim_thread.start()

	main_view = ui.load_view("ising_2d_ios")

	scene_view = scene.SceneView()
	scene_view.frame = (16, 90, 512, 512)
	lattice_scene = LatticeScene(2**int(main_view['lattice_size_slider'].value) * 4 - 1)
	scene_view.scene = lattice_scene
	main_view.add_subview(scene_view)

	plotter_view = scene.SceneView()
	plotter_view.frame = (16, 90, 480, 512)
	plotter_scene = PlotterScene()
	plotter_view.scene = plotter_scene
	main_view['running_view'].add_subview(plotter_view)

	main_view.present('fullscreen')
	on_lattice_size_changed(main_view['lattice_size_slider'])

	while True:
		pass