import math, random

De = 0.144
Re = 4.4616
alpha = 0.49446
k = 0.00000316681518

PARTICLES = 13
INITIAL_STEP_SIZE = 0.03125
TEMPERATURE = 50
STEPS = 100000

def quadratic_potential(x):
	return 0.5 * x**2

def morse_potential(x):
	return De * math.expm1(-alpha * (x - Re))**2

"""
	new_potentials will be in the form
	[
		p1,
		p2,
		...,
		p(index_moved-1),
		0,
		p(index_moved+1),
		...,
		pn
	]
	with the first index always less than the second

	positions is expected to be in the form
	(
		[x1, y1, z1],
		[x2, y2, z2],
		[x3, y3, z3],
		...,
		[xn, yx, zn]
	)
"""
def potentials(potential, positions, index_moved):
	new_potentials = []
	for i in xrange(index_moved):
			new_potentials.append(potential(distance(positions[i], positions[index_moved])))
	new_potentials.append(0)
	for i in xrange(index_moved + 1, PARTICLES):
			new_potentials.append(potential(distance(positions[i], positions[index_moved])))
	return new_potentials

def distance(r1, r2):
	assert len(r1) == len(r2)
	return sum((r1[i] - r2[i])**2 for i in range(len(r1)))**0.5

def pick_new_lmax(previous_selections):
	highs = filter(lambda x : x[1] > 0.5, previous_selections.iteritems())
	lows = filter(lambda x : x[1] < 0.5, previous_selections.iteritems())
	if len(highs) == 0:
		return min(previous_selections.keys()) / 2.0
	elif len(lows) == 0:
		return max(previous_selections.keys()) * 2.0
	else:
		return math.sqrt(filter(lambda x: x[1] == max(map(lambda x: x[1], lows)), lows)[0][0] * filter(lambda x: x[1] == min(map(lambda x: x[1], highs)), highs)[0][0])

def walk(steps, potential_function, lmax, file_name = ""):

	number_accepted = 0
	sum_U = 0

	positions = list([0, 0, 0] for i in xrange(PARTICLES))

	# distances = {}
	# for i in xrange(len(positions) - 1):
	# 	for j in xrange(i+1, len(positions)):
	# 		distances[(i, j)] = distance(positions[i], positions[j])

	potential_dict = {}
	for i in xrange(PARTICLES - 1):
		for j in xrange(i+1, PARTICLES):
			potential_dict[(i, j)] = potential_function(distance(positions[i], positions[j]))
	U = sum(potential_dict.values())

	if file_name:
		data_file = open(file_name, "w")

	for step in xrange(steps):
		print "%d%% done\r"%(100.0*step/steps),
		particle_to_move = random.randint(0, len(positions) - 1)
		rx = (2 * random.random() - 1)
		ry = (2 * random.random() - 1)
		rz = (2 * random.random() - 1)
		d = random.random() * lmax

		dx = rx * d / distance([0, 0, 0], [rx, ry, rz])
		dy = ry * d / distance([0, 0, 0], [rx, ry, rz])
		dz = rz * d / distance([0, 0, 0], [rx, ry, rz])

		positions[particle_to_move][0] += dx
		positions[particle_to_move][1] += dy
		positions[particle_to_move][2] += dz

		new_potentials = potentials(potential_function, positions, particle_to_move)
		new_U = U
		for i in xrange(particle_to_move):
			new_U += new_potentials[i] - potential_dict[(i, particle_to_move)]
		for i in xrange(particle_to_move+1, PARTICLES):
			new_U += new_potentials[i] - potential_dict[(particle_to_move, i)]

		if new_U < U or random.random() < math.exp(-(new_U - U) / k / T):
			number_accepted += 1
			for i in xrange(particle_to_move):
				potential_dict[(i, particle_to_move)] = new_potentials[i]
			for i in xrange(particle_to_move+1, PARTICLES):
				potential_dict[(particle_to_move, i)] = new_potentials[i]
			U = new_U

			# for i in xrange(len(positions)):
			# 	if i < particle_to_move:
			# 		distances[(i, particle_to_move)] = distance(positions[i], positions[particle_to_move])
			# 	if i > particle_to_move:
			# 		distances[(particle_to_move, i)] = distance(positions[i], positions[particle_to_move])
		else:
			positions[particle_to_move][0] -= dx
			positions[particle_to_move][1] -= dy
			positions[particle_to_move][2] -= dz

		sum_U += U
		if file_name:
			# data_file.write(str(distances) + "\n") # storing the distances eats up a lot of memory, because data required ~ (number of particles)^2
			data_file.write(str(U) + "\n")
			if step == steps - 1:
				data_file.write(str(positions) + "\n")
			# data_file.write("\n")

	if file_name:
		data_file.close()

	return float(number_accepted) / steps, sum_U / steps

def simulate(steps_per_simulation, potential, file_name = ""):
	previous_lmax_selections = {}
	percent_accepted = 0
	lmax = INITIAL_STEP_SIZE
	while abs(percent_accepted - 0.5) > 0.1:
		if previous_lmax_selections:
			lmax = pick_new_lmax(previous_lmax_selections)
		print "New lmax: %f"%lmax
		percent_accepted, average_U = walk(steps_per_simulation, potential, lmax, file_name)
		previous_lmax_selections[lmax] = percent_accepted
		print "U Average: %f"%average_U

def plot(file_name):
	import matplotlib.pyplot as plt
	from mpl_toolkits.mplot3d import Axes3D
	import ast

	data_file = open(file_name)

	data_U = []
	final_positions = []

	for line in data_file:
		try:
			data_U.append(float(line))
		except:
			final_positions = ast.literal_eval(line)

	fig = plt.figure()
	ax1 = fig.add_subplot(111, projection = '3d')
	ax1.scatter(map(lambda x:x[0], final_positions), map(lambda x:x[1], final_positions), map(lambda x:x[2], final_positions))
	# ax2 = fig.add_subplot(212)
	# ax2.hist(data_U, bins = 100)
	plt.show()


if __name__ == "__main__":
	T = TEMPERATURE
	simulate(STEPS, morse_potential, "monte_carlo_3D.dat")
	plot("monte_carlo_3D.dat")