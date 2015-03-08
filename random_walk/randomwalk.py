def calculate(dimensions, sample_size, steps):
	import numpy.random as random

	def walk(steps, dimensions = 1):
		position = list(0 for i in xrange(dimensions))
		for step in xrange(steps):
			position[random.randint(0, dimensions)] += (-1 if random.random()<0.5 else 1)
		return tuple(position)

	def magnitude(position):
		return sum(map(lambda x:x**2, position))**0.5

	data_file = open("random_walk_data_%dD.txt"%dimensions, 'w')

	for i in xrange(sample_size):
		if i%10 == 0:
			print '\r%d%% finished'%(100.0*i/sample_size),
		data_file.write(str(magnitude(walk(steps, dimensions))) + '\n')

	data_file.close()

def plot(dimensions):
	import matplotlib.pyplot as plt

	data_file = open("random_walk_data_%dD.txt"%dimensions)

	data = []
	for line in data_file:
		data.append(float(line))

	plt.hist(data, bins = 100)
	plt.show()

if __name__ == "__main__":
	dimensions = 10
	sample_size = 10**5
	steps = 10**3

	calculate(dimensions, sample_size, steps)
	plot(dimensions)
