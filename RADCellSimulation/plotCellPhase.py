import matplotlib.pyplot as plt
import numpy as np
import  math

data = np.genfromtxt('cellPhase.csv', delimiter=',',skip_header=0)
data1 = np.genfromtxt('cellPhase_dose_1.csv', delimiter=',',skip_header=0)
def func_cellPhase(data):
    
    print type(data)
    print data.shape
    time = data[:,0]
    G1Num=data[:,1]
    SNum=data[:,2]
    G2Num=data[:,3]
    num_mitosis = data[:,4]
    num_total = data[:,5]
    G1Ratio=G1Num/num_total
    SRatio=SNum/num_total
    G2Ratio=G2Num/num_total
    mitosisRatio= num_mitosis/num_total
    print mitosisRatio
    
    return (time,G1Ratio,SRatio,G2Ratio,mitosisRatio)



plt.hold(True)
#plt.subplot(2,1,1)
plt.plot(func_cellPhase(data)[0],func_cellPhase(data)[1],'k.',label='G1 Phase Cell')
plt.plot(func_cellPhase(data)[0],func_cellPhase(data)[2],'b.',label='S Phase Cell')
plt.plot(func_cellPhase(data)[0],func_cellPhase(data)[3],'g.',label='G2 Phase Cell')   
plt.plot(func_cellPhase(data)[0],func_cellPhase(data)[4],'r.',label='M Phase Cell')

plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)
plt.xlabel('Time, minutes',fontsize=15)
plt.ylabel('Cell Phase Ratio',fontsize=15)
plt.legend(loc='best')

#plt.subplot(2,1,2)
#plt.plot(func_cellPhase(data1)[0],func_cellPhase(data1)[1],'b.',label='G1 Phase Cell')
#plt.plot(func_cellPhase(data1)[0],func_cellPhase(data1)[2],'g.',label='S Phase Cell')
#plt.plot(func_cellPhase(data1)[0],func_cellPhase(data1)[3],'k.',label='G2 Phase Cell')   
#plt.plot(func_cellPhase(data1)[0],func_cellPhase(data1)[4],'r.',label='M Phase Cell')

#plt.tick_params(axis='x', labelsize=15)
#plt.tick_params(axis='y', labelsize=15)
#plt.xlabel('Time, minutes',fontsize=15)
#plt.ylabel('Cell Phase Ratio',fontsize=15)
#plt.legend(loc='best')


plt.show()

