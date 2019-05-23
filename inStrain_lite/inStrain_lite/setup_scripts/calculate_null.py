import numpy as np
from collections import defaultdict

def simulate(N, p=0.001):

	#Assume "A" is the correct nucleotide.
	result = np.random.choice(['A','G','C','T'], size=N, p=[1-p,p/3,p/3,p/3])
	num_G = np.sum(result == 'G')
	num_C = np.sum(result == 'C')
	num_T = np.sum(result == 'T')

	nucleotides = [num_G,num_C,num_T]
	return nucleotides





def main():
	#Set N
	for N in range(5,2000, 1):
		#bootstraps
		err_obs = defaultdict(int)

		for bt in range(1, 100000):
			# Do simulation
			res = simulate(N)
			err_obs[res[0]] += 1
			err_obs[res[1]] += 1
			err_obs[res[2]] += 1

		# print(dict(err_obs))
		print(str(N) + "\t" + str(err_obs[0]) + "\t" +  str(err_obs[1]) + "\t" + str(err_obs[2]) + "\t" + str(err_obs[3]) + "\t" + str(err_obs[4]) + "\t" + str(err_obs[5]) + "\t" + str(err_obs[6]) + "\t" + str(err_obs[7]) + "\t" + str(err_obs[8]) + "\t" + str(err_obs[9]) + "\t" + str(err_obs[10]) + "\t" + str(err_obs[11]) + "\t" + str(err_obs[12]) + "\t" + str(err_obs[13]) + "\t" + str(err_obs[14]) + "\t" + str(err_obs[15]) + "\t" + str(err_obs[16]) + "\t" + str(err_obs[17]) + "\t" + str(err_obs[18]) + "\t" + str(err_obs[19]) + "\t" + str(err_obs[20]))
	print("done")
if __name__ == '__main__':
	main()
