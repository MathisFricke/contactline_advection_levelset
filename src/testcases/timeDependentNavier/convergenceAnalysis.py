    # -*- coding: utf-8 -*-

import sys
import os
import numpy as np
from matplotlib import pyplot as plt

#testcases are in the same folder as the script
folder = "./"

DeltaX = []
maxDifferenceAngle = []
maxDifferencePosition = []

theoreticalPlotted = True

#def reference(t):
 #   c1 = 0.1
  #  c2 = -2
   # return np.pi/2 + np.arctan(-1/np.tan(theta0) * np.exp(2*c1*t) + c2 * (np.exp(2*c1*t) -1)/(2*c1))

#def reference(t):
 #   return np.pi/2 + np.arctan(-1/np.sqrt(3) * np.exp(0.2*t) + 10*(np.exp(0.2*t) - 1))
for case in os.listdir(folder):
    if (os.path.isfile(folder + case)):
        continue

    angleData = np.genfromtxt(folder + case + "/contactAngle.csv", delimiter = ",")
    positionData = np.genfromtxt(folder + case + "/position.csv", delimiter = ",")
    
    
    #if (not theoreticalPlotted):
     #   plt.plot(curvatureData[:, 0], curvatureData[:, 2])
      #  plt.xlabel("Time in arbitrary units")
        #plt.ylabel("Curvature")
       # theoreticalPlotted = True
    


    DeltaX.append(0.1/float(case[7:]))
    
    timesteps = []
    theoreticalAngle = []
    actualAngle = []

    theoreticalPosition = []
    actualPosition = []
    
    try:
        length = len(angleData[:, 0])
    except IndexError:
        length = 1
        
    for i in range(length):
        time = angleData[i, 0]
        timesteps.append(time)
        actualAngle.append(angleData[i, 1])
        theoreticalAngle.append(angleData[i, 2])
        actualPosition.append(positionData[i, 1])
        theoreticalPosition.append(positionData[i, 2]) 	

    theoreticalAngle = np.array(theoreticalAngle)
    actualAngle = np.array(actualAngle)
    
    theoreticalPosition = np.array(theoreticalPosition)
    actualPosition = np.array(actualPosition)

    maxDiffAngle = max(np.abs(theoreticalAngle - actualAngle))
    maxDiffPosition = max(np.abs(theoreticalPosition - actualPosition))
    
    maxDifferenceAngle.append(maxDiffAngle)
    maxDifferencePosition.append(maxDiffPosition)

print(DeltaX)
print(maxDifferenceAngle)
print(maxDifferencePosition)

plt.subplot(2, 1, 1)
plt.plot(DeltaX, maxDifferenceAngle, linestyle="None", marker ="x")
plt.xlabel("Cell width")
plt.ylabel("max difference in angle")
#plt.savefig("maxDiffAngle.pdf")

plt.subplot(2, 1, 2)
plt.plot(DeltaX, maxDifferencePosition, linestyle="None", marker="x")
plt.xlabel("Cell width")
plt.ylabel("max differnce in position")
plt.show()     
