import matplotlib.pyplot as plt
import numpy as np


data = np.genfromtxt('testGaussian.csv', delimiter=',',skip_header=0)

print type(data)
print data.shape
gaussian_numbers = data[:,0]
print gaussian_numbers
plt.title("Gaussian Histogram")
plt.xlabel("Value")
plt.ylabel("Frequency")

n, bins, patches = plt.hist(gaussian_numbers, 100, normed=1, facecolor='green', alpha=0.75)

plt.show()
