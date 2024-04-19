import math
import numpy as np
import matplotlib.pyplot as plt
plt.rc('font', family='serif')
import os
import glob

folderPath = glob.glob('cav_probes_line1/*')
latestTime = max(folderPath, key=os.path.getctime)
latestTime = os.path.basename(latestTime)

# -------------- Load Data --------------------#
T0 = np.loadtxt(os.path.join('cav_probes_line1/',latestTime,'line_T_p.xy'));
T45 = np.loadtxt(os.path.join('cav_probes_line2/',latestTime,'line_T_p.xy'));
T90 = np.loadtxt(os.path.join('cav_probes_line3/',latestTime,'line_T_p.xy'));
T135 = np.loadtxt(os.path.join('cav_probes_line4/',latestTime,'line_T_p.xy'));
T180 = np.loadtxt(os.path.join('cav_probes_line5/',latestTime,'line_T_p.xy'));
T225 = np.loadtxt(os.path.join('cav_probes_line6/',latestTime,'line_T_p.xy'));
T270 = np.loadtxt(os.path.join('cav_probes_line7/',latestTime,'line_T_p.xy'));
T315 = np.loadtxt(os.path.join('cav_probes_line8/',latestTime,'line_T_p.xy'));
pTime = np.loadtxt('probes_cav0deg/0/p');

#dataFullT = [T0[:,1],T45[:,1],T90[:,1],T135[:,1],T180[:,1],T225[:,1],T270[:,1],T315[:,1]]
dataFullT = [T0[:,1],T90[:,1],T180[:,1],T270[:,1]]
#dataFullT = [T0[:,1],T90[:,1],T135[:,1],T180[:,1],T225[:,1],T270[:,1],T315[:,1]]

#print(len(T0))
#print(len(T45))
#print(len(T90))
#print(len(T135))
#print(len(T180))
#print(len(T225))
#print(len(T270))
#print(len(T315))

dataAvgT = np.mean(dataFullT, axis=0)
rArray = T0[:,0]+0.125

# -------------- Plot Data --------------------#
plt.figure(1)
plt.plot(rArray,dataAvgT,'k-',label='Core Average')
plt.plot((T0[:,0]+0.125),T0[:,1],'r-',label='0deg Instantaneous')
plt.xlabel('Temperature (K)')
plt.ylabel('Radius (m)')
plt.legend(loc="lower right")
#plt.xlim([0 1])
#plt.ylim([0 1])
plt.grid()
plt.savefig('coreT.png')

plt.figure(2)
plt.plot(pTime[:,0],pTime[:,7],'k-',label='Transient Pressure')
plt.xlabel('Time (s)')
plt.ylabel('Pressure (Pa)')
plt.legend(loc="lower right")
#plt.xlim([0 2])
#plt.ylim([0 1])
plt.grid()
plt.savefig('pTime.png')

