#coding:utf8
f = "Input/kmers.txt"
L = []
h = open(f, "rU")
for line in h:
	L.append(line.strip())

print L
L.sort()
print L

for l in L:
	print l