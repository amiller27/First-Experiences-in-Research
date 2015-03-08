"""
	Aaron Miller
	26 Feb 2015

	Set of functions for reading in works of literature, analyzing the frequency of letter
	combinations in those works, and producing writing by a random generation based on
	the collected frequencies

	The data is stored in the following files:
		primary_frequencies.dat
		secondary_frequencies.dat
		tertiary_frequencies.dat
"""

import random

class Writer(object):
	"""
		Allows for collection of data and generation of writing from that data
	"""

	def __init__(self):
		super(Writer, self).__init__()
		
		# load databases into memory
		try:
			pfile = open("primary_frequencies.dat", "U")
			self.primary = dict(filter(None, [None if len(x) < 2 else (x[0], int(x[1])) for x in [l.rsplit(" ", 1) for l in pfile.readlines()]]))
			pfile.close()
		except IOError, e:
			self.primary = {}
		try:
			sfile = open("secondary_frequencies.dat", "U")
			self.secondary = dict(filter(None, [None if len(x) < 2 else (x[0], int(x[1])) for x in [l.rsplit(" ", 1) for l in sfile.readlines()]]))
			sfile.close()
		except IOError, e:
			self.secondary = {}
		try:
			tfile = open("tertiary_frequencies.dat", "U")
			self.tertiary = dict(filter(None, [None if len(x) < 2 else (x[0], int(x[1])) for x in [l.rsplit(" ", 1) for l in tfile.readlines()]]))
			tfile.close()
		except IOError, e:
			self.tertiary = {}

	def _add_combination(self, letters, correlation):
		"""
			Adds the combination of letters to the database indicated by correlation

			letters: the sequence of letters of which to increment the frequency
			correlation: either 1, 2, or 3, depending on which database to add to
		"""

		if correlation == 1:
			if letters in self.primary:
				self.primary[letters] += 1
			else:
				self.primary[letters] = 1
		elif correlation == 2:
			if letters in self.secondary:
				self.secondary[letters] += 1
			else:
				self.secondary[letters] = 1
		elif correlation == 3:
			if letters in self.tertiary:
				self.tertiary[letters] += 1
			else:
				self.tertiary[letters] = 1
		else:
			raise ValueError("Invalid correlation value")

	def close(self):
		"""
			Stores databases to file
		"""

		pfile = open("primary_frequencies.dat", 'w')
		pfile.write("".join(x[0] + " " + str(x[1]) + "\n" for x in self.primary.items()))
		pfile.close()
		sfile = open("secondary_frequencies.dat", 'w')
		sfile.write("".join(x[0] + " " + str(x[1]) + "\n" for x in self.secondary.items()))
		sfile.close()
		tfile = open("tertiary_frequencies.dat", 'w')
		tfile.write("".join(x[0] + " " + str(x[1]) + "\n" for x in self.tertiary.items()))
		tfile.close()


	def train(self, filename):
		"""
			Trains the statistical database with data from a work

			filename: the name of the file to read data from
		"""
		
		try:
			input_file = open(filename, "U")
		except IOError, e:
			raise e

		last = ""
		before_last = ""
		while True:
			c = input_file.read(1)
			if c == "":
				break

			if c == "\n":
				c = " "

			if not last:
				self._add_combination(c, 1)
				last = c
			elif not before_last:
				self._add_combination(c, 1)
				self._add_combination(last + c, 2)
				before_last = last
				last = c
			else:
				self._add_combination(before_last + last + c, 3)
				self._add_combination(last + c, 2)
				self._add_combination(last, 1)
				before_last = last
				last = c

		input_file.close()

	def generate_probabilities(self):
		sample_size = sum(self.primary.values())
		self.primary_probs = dict((x[0], float(x[1]) / sample_size) for x in self.primary.iteritems())
		probs = self.primary_probs.items()
		self.primary_distribution = [probs[0]]
		for x in probs[1:]:
			self.primary_distribution.append((x[0], self.primary_distribution[-1][1] + x[1]))

		self.secondary_probs = {}
		self.secondary_distribution = {}
		for start_string in set(map(lambda x: x[:-1], self.secondary.keys())):

			sample_size = sum(x[1] for x in filter(lambda x: x[0][:-1] == start_string, self.secondary.iteritems()))

			self.secondary_probs[start_string] = dict((x[0][-1], float(x[1]) / sample_size) for x in filter(lambda x: x[0][:-1] == start_string, self.secondary.iteritems()))

			probs = self.secondary_probs[start_string].items()
			
			self.secondary_distribution[start_string] = [(probs[0][0][-1], probs[0][1])]

			for x in probs[1:]:
				self.secondary_distribution[start_string].append((x[0][-1], self.secondary_distribution[start_string][-1][1] + x[1]))

		self.tertiary_probs = {}
		self.tertiary_distribution = {}
		for start_string in set(map(lambda x: x[:-1], self.tertiary.keys())):

			sample_size = sum(x[1] for x in filter(lambda x: x[0][:-1] == start_string, self.tertiary.iteritems()))

			self.tertiary_probs[start_string] = dict((x[0][-1], float(x[1]) / sample_size) for x in filter(lambda x: x[0][:-1] == start_string, self.tertiary.iteritems()))

			probs = self.tertiary_probs[start_string].items()

			self.tertiary_distribution[start_string] = [(probs[0][0][-1], probs[0][1])]

			for x in probs[1:]:
				self.tertiary_distribution[start_string].append((x[0][-1], self.tertiary_distribution[start_string][-1][1] + x[1]))

	def write(self, length):
		result = ""
		# print "{"
		# for i in self.primary_probs.iteritems():
		# 	print "\t" + str(i)
		# print "}"
		# print
		# print "{"
		# for i in self.secondary_probs.iteritems():
		# 	print "\t" + str(i[0]) + " {"
		# 	for j in i[1].iteritems():
		# 		print "\t\t" + str(j)
		# 	print "\t" + "}"
		# print "}"
		# print
		# for i in self.tertiary_probs.iteritems():
		# 	print "\t" + str(i[0]) + " {"
		# 	for j in i[1].iteritems():
		# 		print "\t\t" + str(j)
		# 	print "\t" + "}"
		# print "}"
		# print

		# print "["
		# for i in self.primary_distribution:
		# 	print "\t" + str(i)
		# print "]"
		# print
		# print "{"
		# for i in self.secondary_distribution.iteritems():
		# 	print "\t" + str(i[0]) + " ["
		# 	for j in i[1]:
		# 		print "\t\t" + str(j)
		# 	print "\t" + "]"
		# print "}"
		# print
		# print "{"
		# for i in self.tertiary_distribution.iteritems():
		# 	print "\t" + str(i[0]) + " ["
		# 	for j in i[1]:
		# 		print "\t\t" + str(j)
		# 	print "\t" + "]"
		# print "}"
		# print

		for i in range(length):

			r = random.random()
			
			if len(result) > 1 and result[-2:] in self.tertiary_distribution:
				for i in range(len(self.tertiary_distribution[result[-2:]])):
					if self.tertiary_distribution[result[-2:]][i][1] > r:
						result += self.tertiary_distribution[result[-2:]][i][0]
						break

			elif len(result) > 0 and result[-1] in self.secondary_distribution:
				for i in range(len(self.secondary_distribution[result[-1]])):
					if self.secondary_distribution[result[-1]][i][1] > r:
						result += self.secondary_distribution[result[-1]][i][0]
						break

			else:
				for i in range(len(self.primary_distribution)):
					if self.primary_distribution[i][1] > r:
						result += self.primary_distribution[i][0]
						break

		return result



if __name__ == "__main__":
	from sys import argv
	w = Writer()
	w.train(argv[1])
	w.generate_probabilities()
	print w.write(1000)
	w.close()