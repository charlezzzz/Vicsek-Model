import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import moment, kurtosis, tstd
import math
from math import sqrt, pow



'''	Binder Cumulant:

		The kurtosis is the fourth standardized moment, defined as
				Kurt[X] = E[(X-mu)^4] / (E[(X-mu)^2])^2
		This moment is a measure of the "tailedness" of a distribution 
		as values ofXcloseto the mean, once raised to the fourth power, 
		will add negligible contributions to Kurt[X]. 

		For a Gaussian distribution 
			N(x|mu, sigma^2) = 1/sqrt[2pi*sigma] * e^[-(x-mu)^2/sigma^2]
		the fourth moment is 3*sigma^4, and so Kurt[N]=3. 

		The "Kurtosis excess", g, is then defined using the Gaussian distribution
		as a reference point, 
					g[X] = Kurt[X] - 3
		such that g[N] = 0.

		The "Binder Cumulant", G, is similarly defined as
					G[X] = -g[X]/3 = 1 - Kurt[X]/3
		In the disordered-phase of the Vicsek model, veocity components are random 
		vectors; if drawn from a zero-mean Gaussian, one would expect G[phi]=1/3 in 
		a 2D system. In the ordered-phase, it is expected that Kurt[phi]=sigma[phi]^4
		and so G[phi]=2/3
'''
class CalculateBinderCumulant:

	def __init__(self, data=None):
		self.data = data
		print("\n\nCalculateBinderCumulant object created:\n")


	def make_distribution(self, bins):
		np_array = np.asarray(self.data)
		hist, bin_edges = np.histogram(np_array, bins=bins, normed=True)
		distribution = self.normalize_distribution(hist)
		return distribution, bin_edges


	''' Method:     "normalize_distribution"
		input(s):    A numpy list (p_x) that is a non-normalized distribution
		output(s):   A numpy list (n_x) that is normalized
		description: The normalization condition for a distribution p(x_i) is
							sum_i p(x_i) = 1
					 This method calculates
					 		sum_i p(x_i) = A
					 and then normalizes the distribution through
					 		n(x_i) = p(x_i) / A
					 such that
					 		sum_i n(x_i) = 1
	'''
	def normalize_distribution(self, p_x):
		A = p_x.sum()
		n_x = [None]*len(p_x)
		for i in range(0, len(p_x)):
			n_x[i] = p_x[i] / A
		return np.asarray(n_x)

	def calculate_binder_cumulant(self, distribution):
		kurt = kurtosis(a=distribution, fisher=False)
		binder_cumulant = 1 - kurt/3
		return binder_cumulant

	def test_binder_cumulant(self, distribution):
		kurt = kurtosis(a=distribution, fisher=False)
		binder_cumulant = 1 - kurt/3
		return binder_cumulant


	def gaussian_driver(self):
		# loc=mean, scale=standard deviation, size=num data points
		gausssian = np.random.normal(loc=1.00, scale=2.0, size=1000)

		sns.set(color_codes=True)
		sns.distplot(gausssian)
		#plt.show()

		std_dev = tstd(a=gausssian)
		variance = moment(a=gausssian, moment=2)
		kurt = kurtosis(a=gausssian)
		kurt_fish = kurtosis(a=gausssian, fisher=False)
		print("sigma = %f" %(std_dev))
		print("sigma^2 = %f" %(variance))
		print("Kurt = %f" %(kurt))
		print("kurt_f = %f" %(kurt_fish))


def main():
	calc_bc = CalculateBinderCumulant()
	calc_bc.gaussian_driver()
	


if __name__ == "__main__":
	main()
