f = open('metadata.txt')
samples = {}
for line in f.readlines():
	sample = line.split()[0]
	samples[sample] = ",".join(line.split("\t"))
f.close()

f = open('Cdb.csv')

print(f.readline().strip() + ",sample,Treat_Control,Time_Point,Plot,Depth")
for line in f.readlines():
	s = "_".join(line.split(",")[0].split("_")[:4])
	print(line.strip() + "," + samples[s])
