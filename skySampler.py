#KinMS Sample clouds from sky distribution
"""
Copyright (C) 2017, Mark D. Smith
E-mail: mark.smith -at- physics.ox.ac.uk

Updated versions of the software are available through github:
https://github.com/Mark-D-Smith/KinMS-skySampler

This is a plugin for KinMS, please note the following from that code:
    If you have found this software useful for your research,
    I would appreciate an acknowledgment to the use of the
    "KINematic Molecular Simulation (KinMS) routines of Davis et al., (2013)".
    [MNRAS, Volume 429, Issue 1, p.534-555]

This software is provided as is without any warranty whatsoever.
For details of permissions granted please see LICENCE.md
"""


import numpy as np


def skySampler(sb, inc, posAng, phaseCentre=np.array([0,0,0]),cloudResolution = 20, maxCloudComponents = 0, cloudCap = 3000000, cellSize = 1., sbMode = 'list', clipLevel = 0.):
    """
    
    A function to produce clouds modelling an arbitrary sky distribution. 
    Takes inputs specifing the distribution and observed parameters.
    Returns the list of clouds taken from this model.

    Parameters
    ----------
    
    sb : either list [x(px),y(px),I], Moment 0 map [[I,...,I],[],...,[]] or full datacube (e.g. a CASA .model file) - must be made of un-restored CLEAN components
    
    inc : float
        The kinematic inclination, currently must be single valued
        
    posAng : float
        The Kinematic Position Angle, currently must be single valued
        
    phaseCentre : 2-vector describing the morphological centre of the profile relative to the centre pixel in units of px

    Other Parameters
    ----------------
    
    sbMode : str, optional
        Allows the user to select if a list ('list') of CLEAN components [x,y,I], a Moment 0 map ('mom0') or 
        an observational cube ('cube') has been provided
    
    cloudResolution: int, optional
        Allows the user to specify the number of cloudlets that will be produced for the dimmest spaxel in the cube
        computation time will scale with this value
    
    maxCloudComponents: int, optional
        Sets the maximum number of clouds that should be generated for a particularly bright component,
        if set to 0 (as default), does not limit the maximum number of clouds
        
    cloudCap: int, optional
        Sets the maximum number of clouds that the sampler is allowed to use
    
    cellSize : float, optional
        Sets the number of arcseconds per pixel in making the conversion from pixel units to arcseconds
        
    clipLevel : float, optional
        Set to require the program to clip the noise at the specified level.
    
    Returns
    -------
    
    clouds : the list of cloud positions in the galaxy plane in the format [x("),y("),z("),I]
        the KinMS inputs are then inClouds = clouds[:,0:3] and flux_clouds=clouds[:,3] 
    
    """
    
    
    #Step 1: Convert Cube ==> Moment 0 ==> List
    if sbMode == 'cube':
        #Project to get a moment 0
        sb=sb.sum(axis=2)
    
    if ((sbMode == 'cube') | (sbMode == 'mom0')):
        #Convert into a list of clouds
        cent = [sb.shape[0]/2 + phaseCentre[0]/cellSize,sb.shape[1]/2 + phaseCentre[1]/cellSize]
        sbList = np.zeros([1,3])
        for i in range(0,sb.shape[0]):
            for j in range(0,sb.shape[1]):
                newCloud = np.array([[i,j,sb[i,j]]])
                
                sbList = np.append(sbList,newCloud,axis=0)
        sb=sbList[1:,:]
    #Clip the noise at the specified level
    if clipLevel:
        keepCount = np.sum(i > clipLevel for i in sb[:,2])
        sbKept = np.zeros([keepCount,3])
        j=0
        for i in range(0,sb.shape[0]):
            if sb[i,2] > clipLevel:
                sbKept[j,:] = sb[i,:]
                j = j + 1
    sb=sbKept
    
    
    inFlux = sb[:,2].sum()
    #Now we know we have a list of clouds [x,y,I] in the sky plane
    
    #Step 2: Slice each list element into multiple cloudlets to allow for velDisp and diskThick
    
    
    #Find minimum non-zero value for intensity - divide this by the cloudResolution to get the approximate standard intensity of a single new cloud
    min = np.nanmin(np.where(sb[:,2]>0,sb[:,2],float('Inf')))
    cloudI = min/cloudResolution
    
    
    # Create a new list of clouds that have this intensity
    sbNew = np.zeros([1,3])
    #Work out how many clouds to use for each component
    clouds2Use = np.where(maxCloudComponents * cloudI > sb[:,2], 
        np.full(sb.shape[0],maxCloudComponents),
        np.rint(sb[:,2]/cloudI))
    
    #Limit the maximum number of clouds to cloudCap
    if np.sum(clouds2Use) > cloudCap:
        clouds2Use = np.floor( clouds2Use * cloudCap / np.sum(clouds2Use))
    
    #Initialise an array of length to contain all these clouds
    clouds = np.zeros([int(np.sum(clouds2Use)),4])
    j=0
    sbNew = np.zeros([int(sb.shape[0])])
    for i in range(0,clouds2Use.shape[0]):
        cloudCounter = clouds2Use[i]
        while cloudCounter > 0:
            newCloud=np.array([[(sb[i,0]-cent[0])*cellSize,(sb[i,1]-cent[1])*cellSize,0*cellSize,sb[i,2]/clouds2Use[i]]])
            clouds[j,:]=newCloud
            j=j+1
            cloudCounter = cloudCounter - 1
    outFlux = clouds[:,3].sum()
    
    if (inFlux-outFlux > (inFlux*1e-8)):
        print('WARNING: Flux may not have been conserved',inFlux,outFlux)
    
    retClouds = clouds
    #Step 5: Output list of new values
    return retClouds