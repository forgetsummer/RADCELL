import matplotlib.pyplot as plt
import numpy as np
import  math

def plotSFCurve(fileName,color,plotName): 
    data = np.genfromtxt(fileName, delimiter=',',skip_header=0)
    print type(data)
    print data.shape
    time = data[:,0]
    totalCellNum=data[:,1]  
    #plt.plot(time,totalCellNum,color+'.-',markersize=10, linewidth=2.0, label=fileName[0:-4])
    plt.plot(time,totalCellNum,color+'.-',markersize=10, linewidth=5.0, label=plotName)
    #plt.yscale('log')
    
    
def plotSFCurve_partTimeRange(fileName,color,plotName): 
    data = np.genfromtxt(fileName, delimiter=',',skip_header=0)
    print type(data)
    print data.shape
    time = data[:,0]
    totalCellNum=data[:,1]  
    time_starting = time[11000:,]
    totalCellNum_starting = totalCellNum[11000:,]
    #plt.plot(time,totalCellNum,color+'.-',markersize=10, linewidth=2.0, label=fileName[0:-4])
    plt.plot(time_starting,totalCellNum_starting,color+'.-',markersize=10, linewidth=5.0, label=plotName)
    #plt.yscale('log')
    
def plotSFCurve_test(fileName,color): 
    data = np.genfromtxt(fileName, delimiter=',',skip_header=0)
    print type(data)
    print data.shape
    time = data[:,0]
    totalCellNum=data[:,1]   
    time_starting = time[11999:,]
    totalCellNum_starting = totalCellNum[11999:,]
    plt.plot(time_starting,totalCellNum_starting,color+'.-',markersize=10, linewidth=5.0, label=fileName[0:-4])
    #plt.plot(time,totalCellNum,color+'.-',markersize=10, linewidth=2.0, label=plotName)
    #plt.yscale('log')
#plt.yticks(np.arange(0,1.1,0.05))  
plt.rcParams['figure.figsize'] = (9.33, 7) # set the default plot size

plt.hold(True)
#plotSFCurve('tumorCellNumber_withoutRadiation.csv','r')
plotSFCurve('cellNumber_multipleRadiation_dose=0Gy.csv','r','Tumor Cell Number with 0 Gy dose')
plotSFCurve('cellNumber_multipleRadiation_dose=1Gy.csv','b','Tumor Cell Number with 1 Gy dose')
plotSFCurve('cellNumber_multipleRadiation_dose=2Gy.csv','g','Tumor Cell Number with 2 Gy dose')
plotSFCurve('cellNumber_multipleRadiation_dose=3Gy.csv','m','Tumor Cell Number with 3 Gy dose')
plotSFCurve('cellNumber_multipleRadiation_dose=4Gy.csv','y','Tumor Cell Number with 4 Gy dose')
plotSFCurve('cellNumber_multipleRadiation_dose=5Gy.csv','c','Tumor Cell Number with 5 Gy dose')
plotSFCurve('cellNumber_multipleRadiation_dose=6Gy.csv','k','Tumor Cell Number with 6 Gy dose')
#plotSFCurve('fractionDose_MRT_10Gy_fiveDose_Gap_1000MCS_cellNumber.csv','y','Tumor Cell Number with 10 Gy dose')

plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)
plt.xlabel('Time, MCS',fontsize=15)
plt.ylabel('Tumor Cell Number',fontsize=15)
plt.grid(True,which="both",ls="-")
plt.legend(loc='best')
plt.savefig('tumorCellNumberVSTime.svg')


##############################################################################################
plt.figure()
plotSFCurve('single_MRT_8Gy_5fractions_cellNumber.csv','b','40Gy in 5 fractions')
plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)
plt.xlabel('Time, MCS',fontsize=15)
plt.ylabel('Tumor Cell Number',fontsize=15)
plt.grid(True,which="both",ls="-")
plt.legend(loc='best')
plt.savefig('tumorCellNumberVSTime_fractionateDose_8Gy.svg')

plt.figure()
plotSFCurve('single_MRT_20Gy_2fractions_cellNumber.csv','g','hyper-fractionated: 40Gy in 2 fractions ')
plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)
plt.xlabel('Time, MCS',fontsize=15)
plt.ylabel('Tumor Cell Number',fontsize=15)
plt.grid(True,which="both",ls="-")
plt.legend(loc='best')
plt.savefig('tumorCellNumberVSTime_fractionateDose_10Gy.svg')

plt.figure()
plotSFCurve_partTimeRange('single_MRT_8Gy_5fractions_cellNumber.csv','b','40Gy in 5 fractions')
plotSFCurve_partTimeRange('single_MRT_20Gy_2fractions_cellNumber.csv','g', '40Gy in 2 fractions ')
plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)
plt.xlabel('Time, MCS',fontsize=15)
plt.ylabel('Tumor Cell Number',fontsize=15)
plt.grid(True,which="both",ls="-")
plt.legend(loc='best')
plt.savefig('tumorCellNumberVSTime_fractionateDose.svg')



#############################################################################################################
plt.figure()
plotSFCurve_partTimeRange('fractionDose_MRT_8Gy_fiveDose_Gap_1000MCS_cellNumber.csv','b','Multi array MRT, 40Gy in 5 fractions')
plotSFCurve_partTimeRange('fractionDose_MRT_10Gy_fiveDose_Gap_1000MCS_cellNumber.csv','m','Multi array MRT,50Gy in 5 fractions')
plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)
plt.xlabel('Time, MCS',fontsize=15)
plt.ylabel('Tumor Cell Number',fontsize=15)
plt.grid(True,which="both",ls="-")
plt.legend(loc='best')
plt.savefig('tumorCellNumberVSTime_MRT_fractionateDose.svg')

plt.figure()
plotSFCurve('fractionDose_MRT_8Gy_fiveDose_Gap_1000MCS_cellNumber.csv','m','Multi array MRT,40Gy in 5 fractions')
plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)
plt.xlabel('Time, MCS',fontsize=15)
plt.ylabel('Tumor Cell Number',fontsize=15)
plt.grid(True,which="both",ls="-")
plt.legend(loc='best')
plt.savefig('tumorCellNumberVSTime_MRT_fractionateDose_8Gy.svg')

plt.figure()
plotSFCurve('fractionDose_MRT_10Gy_fiveDose_Gap_1000MCS_cellNumber.csv','y','Multi array MRT,50Gy in 5 fractions')
plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)
plt.xlabel('Time, MCS',fontsize=15)
plt.ylabel('Tumor Cell Number',fontsize=15)
plt.grid(True,which="both",ls="-")
plt.legend(loc='best')
plt.savefig('tumorCellNumberVSTime_MRT_fractionateDose_8Gy.svg')


plt.show()
