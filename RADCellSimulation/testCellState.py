import matplotlib.pyplot as plt
import numpy as np
import  math

data = np.genfromtxt('cellState.csv', delimiter=',',skip_header=0)

print type(data)
print data.shape
time = data[:,0]
S1Num=data[:,1]
S21Num=data[:,2]
S22Num=data[:,3]
S3Num = data[:,4]
num_total =S1Num + S21Num + S22Num + S3Num
S1Ratio=S1Num/num_total
S21Ratio=S21Num/num_total
S22Ratio=S22Num/num_total
S3Ratio= S3Num/num_total
SF = (S1Num+S21Num)/num_total
print 'the sf is '
print SF
print 'the size of sf is', len(SF)
plt.hold(True)
time = time*0.01
plt.plot(time,S1Ratio,'b.',label='S1 State Cell')
plt.plot(time,S21Ratio,'g.',label='S21 State Cell')
plt.plot(time,S22Ratio,'k.',label='S22 State Cell')   
plt.plot(time,S3Ratio,'r.',label='S3 State Cell')

plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)
plt.xlabel('Dose, Gy',fontsize=15)
plt.ylabel('Cell State Ratio',fontsize=15)
#plt.xlabel('Time, minutes',fontsize=15)
#plt.ylabel('Cell State Ratio',fontsize=15)
plt.legend(loc='best')
plt.figure()
plt.plot(time,SF,'b.',label='Cell Survival Fraction')
#plt.yscale('log')
plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)
plt.xlabel('Dose, Gy',fontsize=15)
plt.ylabel('SF',fontsize=15)
#plt.xlabel('Time,minutes',fontsize=15)
#plt.ylabel('SF',fontsize=15)
plt.legend(loc='best')

def  plotSFCurveWRTCellPhase(data):
    time = data[:,0]
    S1Num=data[:,1]
    S21Num=data[:,2]
    S22Num=data[:,3]
    S3Num = data[:,4]
    num_total =S1Num + S21Num + S22Num + S3Num
    SF = (S1Num+S21Num)/num_total
    
    return  SF
    
#data_G1 = np.genfromtxt('cellState_G1.csv', delimiter=',',skip_header=0)
#data_G2 = np.genfromtxt('cellState_G2.csv', delimiter=',',skip_header=0)
#data_S = np.genfromtxt('cellState_S.csv', delimiter=',',skip_header=0)

#plt.figure()

#plt.plot(time,plotSFCurveWRTCellPhase(data_G1),'b.',label='Cell Survival Fraction, G1')
#plt.yscale('log')
#plt.hold(True)
#plt.plot(time,plotSFCurveWRTCellPhase(data_G2),'r.',label='Cell Survival Fraction, G2')
#plt.plot(time,plotSFCurveWRTCellPhase(data_S),'g.',label='Cell Survival Fraction, S')
#plt.tick_params(axis='x', labelsize=15)
#plt.tick_params(axis='y', labelsize=15)
#plt.xlabel('Dose, cGy',fontsize=15)
#plt.ylabel('SF',fontsize=15)
##plt.xlabel('Time,minutes',fontsize=15)
##plt.ylabel('SF',fontsize=15)
#plt.legend(loc='best')




plt.show()

