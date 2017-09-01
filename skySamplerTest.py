#Test if skySampler actually works!!!
import sys

import numpy as np

from KinMS import *
from skySampler import *
from makeplots import *
import mge_vcirc as vc



distance =16.5

#1. Create a new model cube using a simulated galaxy
obspars={}
obspars['cellsize']=1.  #arcsec/px
obspars['xsize']=64.    #arcsec
obspars['ysize']=64.    #arcsec
obspars['vsize']=1024.    #km/s
obspars['dv']=1.
obspars['beamsize']=np.array([3.,3.,0.])

params={}
params['inc'] = 70.
params['rflat']=24.
params['rtrunc']=3.
params['vPosAng']=0.
params['posAng']=37.
params['vPhaseCen'] = np.array([0.,0.])    #arcsec
params['phaseCen'] = np.array([0.,0.])

#Previously determined M/L
obspars['psiOuter']=8.25
obspars['rOuter']=5.09
obspars['gasSigma'] = 2.2 #km/s

rad=np.arange(0.,3*obspars['xsize'],0.5)
velrad = rad
velrad[0]=velrad[1]/10**30
velprof=rad*0.

param=np.zeros(9)
param[0] = 6.59          #PsiCent
param[1] = 75          #inc
param[2] = 8.17          #logMBH
param[3] = 0.5          #intFlux
param[4] = 1.43          # r_break
param[5] = 93.2          #vPosAng
param[6] = -0.17         #vPhaseCenX
param[7] = 0.04          #vPhaseCenY
param[8] = -5.73         #vOffset


####################Generate a composite velocity curve###################################

mgeComponents = {}
mgeLogTotalLuminosity = np.array([4.325,3.103,3.659,3.455,3.158,2.818,2.286,1.424,0.842])        #log10 Lsol pc^-2
mgeLogWidth = np.array([-0.839,-0.262,-0.125,0.333,0.684,1.129,1.694,2.026,2.026])               #arcsec
mgeComponents['TotalLuminosity'] = 10 ** mgeLogTotalLuminosity
mgeComponents['Width'] = 10 ** mgeLogWidth
mgeComponents['QSky'] = np.array([0.675,0.950,0.400,0.700,0.781,0.615,0.400,0.400,0.950])                      #'axis ratio'
print mgeComponents['TotalLuminosity']
print mgeComponents['Width']

velBH = vc.mge_vcirc(np.array([0.]), np.array([1.]),np.array([1.]),45., 1., distance, velrad)                  # inclination is irrelevant to this calculation

incVector = np.arange(66.6,89., 1.)
obspars['incVector'] = incVector
velStar=np.zeros((incVector.shape[0],velrad.shape[0]))
i=0
print 'Generating velocity curves'
for inc in incVector:
    velStar[i] = vc.mge_vcirc(mgeComponents['TotalLuminosity'],mgeComponents['Width'],mgeComponents['QSky'],inc,0.,distance,velrad)
    i=i+1
velocityCurves=(velBH,velStar)

### Scale JAM's velocity profiles to correct parameter values
velBH = velocityCurves[0] * np.sqrt( 10 ** 8 )
##### Choose correct velocity curve for inclination:
inc = params['inc']
incVector = obspars['incVector']
myVelCurve = np.zeros(velrad.shape)
for i in np.arange(0,velocityCurves[1].shape[1]):
	velCurveInterFunc = interpolate.interp1d(incVector,velocityCurves[1][:,i],kind='linear')
	myVelCurve[i] = velCurveInterFunc(inc)
velStar = myVelCurve
#print velStar


#psi is M/L -- generate a vector of M/L for the distances
psi = np.zeros(velStar.shape)
for i in np.arange(0,velrad.shape[0]):
    r = velrad[i]
    if r<=param[4]:
        psi[i] = param[0]
    elif (param[4] < r < obspars['rOuter']) :
        m2 = (obspars['psiOuter']-param[0])/(obspars['rOuter']-param[4])
        c2 = param[0] - m2 * param[4]
        psi[i] = m2 * r + c2
    else: psi[i] = obspars['psiOuter']

velStar = np.multiply(velStar, np.sqrt(psi))
velprof = np.sqrt( velBH ** 2 + velStar ** 2 )

##########################################################################################




sbprof=rad*0.
sbprof=np.exp(-rad/params['rflat'])
sbprof[np.where(rad<params['rtrunc'])]=0.

cube=KinMS(obspars['xsize'],obspars['ysize'],obspars['vsize'],obspars['cellsize'],obspars['dv'],obspars['beamsize'],params['inc'],sbMode='axisymmetric',sbProf=sbprof,
    sbRad=rad,velRad=rad,velProf=velprof,fixSeed=True,fileName='model',vPhaseCen=params['vPhaseCen'],phaseCen=params['phaseCen'],
    posAng=params['posAng'],vPosAng=params['vPosAng'],nSamps=3000000,cleanOut=True)

makeplots(obspars['xsize'],obspars['ysize'],obspars['vsize'],obspars['cellsize'],obspars['dv'],obspars['beamsize'],posAng=params['posAng'],overcube=False,pvdthick=1)

cube2=KinMS(obspars['xsize'],obspars['ysize'],obspars['vsize'],obspars['cellsize'],obspars['dv'],obspars['beamsize'],params['inc'],sbMode='axisymmetric',sbProf=sbprof,
    sbRad=rad,velRad=rad,velProf=velprof,fixSeed=True,fileName='model2',vPhaseCen=params['vPhaseCen'],phaseCen=params['phaseCen'],
    posAng=params['posAng'],vPosAng=params['vPosAng'],nSamps=3000000,cleanOut=False)

#cube += np.random.normal(0.,1.,size=cube2.shape)


"""
#For Reference, the KinMS parameters are:
KinMS(xs,ys,vs,cellSize,dv,beamSize,inc, sbMode = 'axisymmetric', samplerMode = 'list', 
    gasSigma=0,sbProf=[],sbRad=[],velRad=[],velProf=[],fileName=False,diskThick=0,cleanOut=False,
    ra=0,dec=0,nSamps=100000,posAng=0.0,intFlux=0,inClouds=[],vLOS_clouds=[],flux_clouds=0,vSys=0,
    restFreq=115.271e9,phaseCen=np.array([0.,0.]),vOffset=0,fixSeed=False,vRadial=0,vPosAng=0,vPhaseCen=np.array([0.,0.]),
    returnClouds=False,gasGrav=False)
"""

#2. Pass to skySampler to generate clouds

"""
#For reference, the skySampler parameters are:
skySampler(sb, inc, posAng, phaseCentre=np.array([0,0,0]),cloudResolution = 20, maxCloudComponents = 0, cloudCap = 3000000, cellSize = 1., sbMode = 'list', clipLevel = 0.):
"""
clouds=skySampler(cube,sbMode='cube',cellSize = obspars['cellsize'])


#3. Feed clouds to KinMS to create a model of the galaxy
newCube=KinMS(obspars['xsize'],obspars['ysize'],obspars['vsize'],obspars['cellsize'],obspars['dv'],obspars['beamsize'],params['inc'],sbMode='skyProfile',
    velRad=rad,velProf=velprof,inClouds=clouds[:,0:3],flux_clouds=clouds[:,3],fixSeed=True,fileName='ssCube',vPhaseCen=np.array([0.,0.]),phaseCen=np.array([0.,0.]),
    posAng=0.,vPosAng=37.)
    

#4. Plot Results
"""makeplots(f,xsize,ysize,vsize,cellsize,dv,beamsize,posang=0,overcube=False,pvdthick=5,nconts=11., rms = 1.,**kwargs)"""
makeplots(cube2,obspars['xsize'],obspars['ysize'],obspars['vsize'],obspars['cellsize'],obspars['dv'],obspars['beamsize'],
    posAng=params['posAng'],overcube=newCube)

