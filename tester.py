# Tester.py for skySampler

from skySampler import *
import matplotlib.pyplot as plt
from KinMS import *
import numpy as np
from makeplots import makeplots
import time

#make list of clouds
print '+++Making list of clouds+++'
t0=time.time()
list = np.zeros((101,101))
centre = [50,50]
for i in range(0,list.shape[0]):
    for j in range(0,list.shape[1]):
        rad = np.sqrt((i-centre[0])**2+((j-centre[1])*0.7)**2)
        if ((rad < 15) & (rad > 5)):
            list[i,j] = 200* np.exp(-(rad/10)**2)
t1 = time.time()
print 'Making clouds took: ' + str(t1-t0) + ' seconds.'

#apply skySampler
clouds = skySampler(list,30,60,sbMode='mom0',cellSize = 0.1)
t2 = time.time()
print 'Sampling the sky took: ' + str(t2-t1) + ' seconds.'

print clouds
#Plot resultant distribution
velRad = np.arange(0,100,0.1)
f = KinMS(10,10,31,0.1,0.1,[0.01,0.01,0],30, posAng = 60, sbMode = 'skyProfile', sbProf = list, samplerMode = 'mom0', velRad = velRad,velProf = velRad ** 0.5, cleanOut = True, fileName = 'testCube.fits', gasSigma = 1)

t3 = time.time()
print 'Making KinMS cube took: ' + str(t3-t2) + ' seconds.'

makeplots(f,100,100,101,0.1,0.1,[0.01,0.01,0],posang = 60)

