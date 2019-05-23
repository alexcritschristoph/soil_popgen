import glob

from collections import defaultdict
data = defaultdict(list)
for fn in glob.glob('./null*.txt'):
	f = open(fn)
	for line in f.readlines():
		if 'done' not in line:
			datum = [int(x) for x in line.strip().split()[1:]]
			num = line.strip().split()[0]
			data[num].append(datum)

for d in data:
	print(d + "\t" + "\t".join([str(sum(x)/1000000) for x in zip(*data[d])]))
