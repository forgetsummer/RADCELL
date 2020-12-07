import matplotlib.pyplot as plt
import numpy as np
import  math



data = np.genfromtxt("microdosimetry_Eng_1PNum_10000_dose.csv", delimiter=',',skip_header=1)
cellDose = data[:,1]

print cellDose
n, bins, patches = plt.hist(cellDose, 50, normed=1, facecolor='green', alpha=0.75)

plt.show()
