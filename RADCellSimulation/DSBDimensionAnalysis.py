import numpy as np 
import matplotlib.pyplot as plt
import matplotlib as mpl

import numpy as np
from scipy.optimize import curve_fit
mpl.rc("font", family="Times New Roman"); plt.rcParams.update({'mathtext.default':  'regular' })

data = np.genfromtxt('Co-60-2MeV_100_DSBDimension.csv', delimiter=',',skip_header=1)

Dimension=data[:,1] 

print Dimension

n, bins, patches = plt.hist(Dimension, 50, normed=1, facecolor='green', alpha=0.75)

plt.show()
