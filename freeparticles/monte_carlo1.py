import math, random

De = 0.144
x0 = 4.4616
alpha = 0.49446
k = 0.00000316681518

def quadratic_potential(x):
	return 0.5 * x**2

def particle_potential(x):
	return De * math.expm1(-alpha * (x - x0))**2

def pick_new_lmax(previous_selections):
	highs = filter(lambda x : x[1] > 0.5, previous_selections.iteritems())
	lows = filter(lambda x : x[1] < 0.5, previous_selections.iteritems())
	if len(highs) == 0:
		return min(previous_selections.keys()) / 2.0
	elif len(lows) == 0:
		return max(previous_selections.keys()) * 2.0
	else:
		return math.sqrt(filter(lambda x: x[1] == max(map(lambda x: x[1], lows)), lows)[0][0] * filter(lambda x: x[1] == min(map(lambda x: x[1], highs)), highs)[0][0])

def walk(steps, potential, lmax, file_name = ""):

	number_accepted = 0
	sum_x = 0
	sum_U = 0

	x = x0
	U = potential(x0)

	if file_name:
		data_file = open(file_name, "w")

	for i in xrange(steps):
		new_x = x + lmax * (2 * random.random() - 1)
		new_U = potential(new_x)

		if new_U < U or random.random() < math.exp(-(new_U - U) / k / T):
			number_accepted += 1
			x = new_x
			U = new_U

		sum_x += x
		sum_U += U
		if file_name:
			data_file.write(str(x) + " " + str(U) + "\n")

	if file_name:
		data_file.close()

	return float(number_accepted) / steps, sum_x  / steps, sum_U / steps

def simulate(steps_per_simulation, potential, file_name = ""):
	previous_lmax_selections = {}
	percent_accepted = 0
	lmax = 1
	while abs(percent_accepted - 0.5) > 0.1:
		if previous_lmax_selections:
			lmax = pick_new_lmax(previous_lmax_selections)
		print "New lmax: %f"%lmax
		percent_accepted, average_x, average_U = walk(steps_per_simulation, potential, lmax, file_name)
		previous_lmax_selections[lmax] = percent_accepted
		print "X Average: %f"%average_x
		print "U Average: %f"%average_U

def plot(file_name):
	import matplotlib.pyplot as plt

	data_file = open(file_name)

	data_x = []
	data_U = []
	for line in data_file:
		data_x.append(float(line.split()[0]))
		data_U.append(float(line.split()[1]))

	fig, axes = plt.subplots(nrows = 2)
	axes[0].hist(data_x, bins = 100)
	axes[1].hist(data_U, bins = 100)
	plt.show()


if __name__ == "__main__":
	T = 5000
	simulate(1000000, particle_potential, "monte_carlo1.dat")
	plot("monte_carlo1.dat")