import matplotlib.pyplot as plt

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
