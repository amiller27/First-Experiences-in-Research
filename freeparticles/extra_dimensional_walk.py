import itertools
import numpy as np

sigma = 1.0
epsilon = 1.0

PARTICLES = 13
CYCLES = 8000
STEPS_PER_WALK = 50
H = 1

def lennard_jones(rsquared):
	sr6 = sigma / rsquared**3
	return 4 * epsilon * (sr6**2 - sr6)

def lj_gradient(positions, index):
	''' generates the gradient of the LJ potential with respect to a given particle '''
	result = list(0 for i in xrange(len(positions[0])))
	for particle in positions:
		for dimension in xrange(len(result)):
			r2 = dsquared(particle, positions[index])
			sr6 = (sigma**2 / r2)**3
			result[dimension] = 24.0 * epsilon * (sr6 - 2 * sr6**2) * (positions[index][dimension] - particle[dimension]) / r2**0.75
	return np.array(result, ndmin=2).T

def random_direction(n):
	''' generates a normalized vector pointing in a random direction in n-dimensional space '''
	vector = np.array([np.random.normal() for i in xrange(dimensions)], ndmin=2)
	length = sum(itertools.imap(lambda x:x**2, vector))**0.5
	return vector.T / length

def normalize(vector):
	''' Returns a normalized copy of the vector '''
	length = sum(itertools.imap(lambda x:x**2, vector.flat))**0.5
	return vector / length

def dsquared(x1, x2):
	assert x1.shape == x2.shape
	return sum((x1[i][0] - x2[i][0])**2 for i in xrange(x1.ndims[0]))

def potentials(potential, positions, index_moved):
	"""
		Calculates ONLY the updated portions of the two-body potential

		result will be in the form
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

		positions is expected to be in the form
		(
			[x1, y1, z1],
			[x2, y2, z2],
			[x3, y3, z3],
			...,
			[xn, yx, zn]
		)
	"""
	new_potentials = []
	for i in xrange(index_moved):
			new_potentials.append(potential(dsquared(positions[i], positions[index_moved])))
	new_potentials.append(0)
	for i in xrange(index_moved + 1, len(positions)):
			new_potentials.append(potential(dsquared(positions[i], positions[index_moved])))
	return new_potentials

def constant_energy_walk(positions, steps):
	i = 0
	for step in xrange(steps):
		grad = lj_gradient(positions, i)
		V = random_direction(len(positions[0]) - 1) + [0.0]

		if not filter(None, grad[:-1]):
			move_direction = V
		else:
			grad_projected = grad[:-1] + [0]
			rotation_90 = [0.5**0.5] + [0]*(len(positions[0]) - 2) + [0.5**0.5]
			axis = quat_rotate(normalize(grad_projected), rotation_90)

		i += 1

def compress(positions):
	pass

def conjugate_minimize(positions):
	pass

def simulate(particles, cycles, steps, extra_dimensions):
	positions = list( \
		list( \
			np.random.uniform(-3 * sigma, 3 * sigma) for j in range(3) \
		) for i in range(particles) \
	)
	positions += [0 for i in range(extra_dimensions)]

	conjugate_minimize(positions)

	for cycle in xrange(cycles):
		constant_energy_walk(positions, steps)
		new_positions = copy_positions(positions)

		# if the compression is not successful, continue from the previous point in D+H dimensions
		if compress(new_positions):
			positions = new_positions
			conjugate_minimize(positions)

	return positions

if __name__ == "__main__":
	simulate(PARTICLES, CYCLES, STEPS_PER_WALK, H)
