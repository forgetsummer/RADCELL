import matplotlib.pyplot as plt
import numpy as np
import  math

def plotSFCurve(fileName,color): 
    data = np.genfromtxt(fileName, delimiter=',',skip_header=0)
    print type(data)
    print data.shape
    dose = data[:,0]
    totalCellNum=data[:,1]
    SF=data[:,2]     
    plt.plot(dose,SF*100,color+'*-',markersize=10, linewidth=2.0, label=fileName[0:-4])
    #plt.yscale('log')
#plt.yticks(np.arange(0,1.1,0.05))  
plt.rcParams['figure.figsize'] = (9.33, 7) # set the default plot size
plotSFCurve('cellState.csv','b')  
#plotSFCurve('cellState_G1_SF.csv','y')
#plotSFCurve('cellState_S_SF.csv','r')
#plotSFCurve('cellSF_without_bystander_effect.csv','r')
#plotSFCurve('cellState_M_SF.csv','k')
plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)
plt.xlabel('Dose, Gy',fontsize=15)
plt.ylabel('Cell Survival Fraction',fontsize=15)
plt.grid(True)
plt.legend(loc='best')
plt.savefig('SF_dose.svg')
plt.show()
