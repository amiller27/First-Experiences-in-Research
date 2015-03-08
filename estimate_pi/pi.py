import matplotlib.pyplot as plt
import numpy.random as random
sample_size = 10**8
batch_size = 100000

plt.ion()

inside_count = 0

for i in xrange(sample_size / batch_size):
	points = list([random.rand(), random.rand()] for j in range(batch_size))
	# inside_x_points = []
	# inside_y_points = []
	# outside_x_points = []
	# outside_y_points = []

	for point in points:
		if point[0]**2 + point[1]**2 <= 1:
			inside_count += 1
	# 		inside_x_points.append(point[0])
	# 		inside_y_points.append(point[1])
	# 	else:
	# 		outside_x_points.append(point[0])
	# 		outside_y_points.append(point[1])
	# plt.scatter(inside_x_points, inside_y_points, s = 2, c = 'r', linewidths = 0)
	# plt.scatter(outside_x_points, outside_y_points, s = 2, linewidths = 0)
	# plt.draw()
	print "Pi is approximately: %f"%(4.0 * inside_count / ((i+1)*batch_size))

plt.ioff()
plt.show()
