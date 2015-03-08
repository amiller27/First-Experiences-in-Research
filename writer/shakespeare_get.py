import urllib2, HTMLParser
import re

text_file = open("shakespeare.txt", 'w')

class HTMLShakespeareParser(HTMLParser.HTMLParser):
	def handle_starttag(self, tag, attrs):
		if tag == "a":
			link = dict(attrs)["href"]
			if "/index.html" in link:
				works.append(link[:-11])

works = []
html_index_file = urllib2.urlopen("http://shakespeare.mit.edu")
HTMLShakespeareParser().feed(html_index_file.read())
html_index_file.close()

for title in works:
	html_data_file = urllib2.urlopen("http://shakespeare.mit.edu/" + title + "/full.html")
	html_data = html_data_file.read()
	html_data_file.close()

	text_file.write(re.sub("\\s+", " ", " " + re.sub("<[^>]*>", " ", html_data.split("<H3>", 1)[1]) + " "))

text_file.close()