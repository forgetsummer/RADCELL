import matplotlib.pyplot as plt
import numpy as np
import  math

def plotSFCurve(fileName,color): 
    data = np.genfromtxt(fileName, delimiter=',',skip_header=0)
    print type(data)
    print data.shape
    time = data[:,0]
    totalCellNum=data[:,1]   
    plt.plot(time,totalCellNum,color+'*-',markersize=10, linewidth=2.0, label=fileName[0:-4])
    plt.yscale('log')
#plt.yticks(np.arange(0,1.1,0.05))  
plt.rcParams['figure.figsize'] = (9.33, 7) # set the default plot size
plotSFCurve('Mono_Electron_Plane_Eng_1PNum_5000000_cellState.csv','b') 

plt.hold(True)
plotSFCurve('tumorCellNumber_withoutRadiation.csv','r')
plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)
plt.xlabel('Time, MCS',fontsize=15)
plt.ylabel('Tumor Cell Number',fontsize=15)
plt.grid(True,which="both",ls="-")
plt.legend(loc='best')
plt.savefig('tumorCellNumberVSTime.svg')
plt.show()
