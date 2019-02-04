# Second build of SkySampler. This program now includes multiple routines - the main, most-costly,
# sampling algorithm and the lightweight transformation program, which rotates the 
# sky-coordinates of clouds to the galaxy plane

# A typical calling order would be to supply the datacube or moment 0, then:
#    clouds = sampleClouds(cube)
#    clouds = transformClouds(clouds, posAng = posAng, inc = inc, cent = cent)
#    clouds = sampleDisc(clouds, discThick, sbRad = rad)
#    modelCube = KinMS(<Blah>, inClouds = clouds[:,0:3], flux_clouds = clouds[:,3])

# Modifications:
# v2. December 2017, MDS, Oxford
# Testing and correction 25 December 2017, MDS, Horsham
# Updated for Python 3 compatibility 26 December 2017, MDS, Horsham

import __future__
import numpy as np
from scipy import interpolate
import time
from copy import copy

import KinMS

def sampleClouds(sb, cellSize, nSamps = 0, sampFact = 20, weighting = None, allow_undersample = False, verbose = True, debug = False):
    """
    Given a 2d image or 3d cube, will generate coordinates for a point-source model of that image
    
    ================
    Inputs:
    
    sb: np.ndarray of 2 or 3 dimensions containing the astronomical data to be modelled. If 3D, the data will be summed along the third axis.
    
    cellSize: float arcsec/px The scale of each pixel, used to convert cloud positions in pixels to physical unity
    
    nSamps: int The total number of point samples to draw from the distribution, where these are distributed weighted by the intensity.
    
    sampFact: int If nSamps not specified, assume a uniform sampling of clouds, with this many per pixel
               WARNING: For large images, this can lead to excessive slow-downs
    
    weighting: np.ndarray of 2 dimensions of same size as sb (or with same first dimensions). Used to weight the sampling by distributions other than intensity, eg. by velocity dispersion
    
    allow_undersample: bool Default FALSE. Prevents using a small number of clouds to sample a large matrix. If the matrix is sparse, this can be disabled.
    
    ================
    Outputs:
    
    clouds: An array of [x,y,I] values corresponding to particle positions relative to the cube centre
    
    """
    
    
    assert len(sb.shape) == 2 or len(sb.shape) == 3, "The input array must be 2 or 3 dimensions"
    
    if len(sb.shape) == 3: sb = sb.sum(axis=2)
    if not nSamps == 0 and allow_undersample == False: assert nSamps > sb.size, "There are insufficiently many clouds to sample the distribution. If this is a sparse array, use allow_undersample = True"
    
    inFlux = sb.sum()
    
    #Convert to a list of pixel values
    cent = [sb.shape[0]/2, sb.shape[1]/2]
    t0 = time.time()
    coords = np.arange(-cent[0],cent[0]) * cellSize
    delY,delX = np.meshgrid(coords,coords)     #Order for consistency with index order
    
    sbList = np.zeros((delX.flatten().shape[0],3))
    sbList[:,0] = delX.flatten()
    sbList[:,1] = delY.flatten()
    sbList[:,2] = sb.flatten()
    sb = sbList
    t1 = time.time()
    
    #Calculate number of clouds to use per pixel. By default, weight uniformly, other modes are intensity-weighted (sampFact set, nSamps not) or custom weighting
    #Priority: Custom, Intensity, Uniform
    if nSamps: 
        if not np.any(weighting):
            scheme = 'Intensity-weighted'
            weighting = sb[:,2]
        else: 
            scheme = 'custom weighting'
            weighting = weighting.flatten()
        weighting=np.abs(weighting)
        intWeight = weighting.sum()
        iCloud = intWeight/nSamps
        nClouds = np.floor(weighting/iCloud)
        
    else: 
        scheme = 'uniform'
        nClouds = np.full(sb[:,2].shape,sampFact)
    if verbose: print('Using a ',scheme,' scheme to sample with ',nClouds.sum(),' clouds.')
    
    # Generate the final list of all clouds
    clouds = np.zeros([int(nClouds.sum()),4])
    k=0
    
    pixSmoothing = 0.5*cellSize              #Smooth cloud positions to within a cell, rather than gridded to the centre of that cell
    for i in np.arange(0,sb.shape[0]):
        if not nClouds[i] == 0:
            for j in np.arange(0,nClouds[i]):
                #print i,j,nClouds[i],k
                clouds[k,:] = np.array([[sb[i,0]+pixSmoothing*np.random.uniform(low=-1.,high=1.),sb[i,1]+pixSmoothing*np.random.uniform(low=-1.,high=1.),0.,sb[i,2]/nClouds[i]]])
                k = k + 1
    t2=time.time()
    if debug:
        print('    Generating pixel positions: ',t1-t0)
        print('    Listing all sub-clouds: ',t2-t1)
    #Sanity checking:
    if not (clouds[:,3].sum() - inFlux) < 1e-3: print('Flux not conserved: '+str(100*(clouds[:,3].sum()-inFlux)/inFlux)+'%')
    #print(clouds)
    return clouds

def transformClouds(inclouds, posAng = 90., inc = 0., cent = [0.,0.]):
    """
    Calculate the galaxy co-ordinates of clouds from the sky plane. This MUST be used if any of the following conditions are true:
    inc != 0
    posAng != 90
    cent != [0,0]
    
    This exists as a stand-alone routine since an MCMC fit to the galaxy will likely need to run this every step, 
    and the sampleClouds routine is computationally expensive
    ============
    Inputs:
    
    clouds: np.ndarray The output of the sampleClouds array [x,y,I]
    
    posAng: 0=<float<360 The position angle of the galaxy major axis, measured from the y-axis of the cube
    
    inc: 0=<float<90 The inclination of the galaxy, with 90 being edge-on
    
    cent: [float,float] The photometric centre of the galaxy relative to the centre of the cube
    
    ============
    Outputs:
    
    clouds: np.ndarray Positions and intensities of particles [x',y',I]
    """
    clouds = copy(inclouds)
    clouds[:,0:2] = np.array([clouds[:,0] - cent[0],clouds[:,1] - cent[1]]).T
    posAng = np.radians(90-posAng)
    xNew = np.cos(posAng) * clouds[:,0] - np.sin(posAng) * clouds[:,1]
    yNew = np.sin(posAng) * clouds[:,0] + np.cos(posAng) * clouds[:,1]
    yNew = yNew / np.cos(np.radians(inc))
    
    clouds[:,0] = xNew
    clouds[:,1] = yNew
    
    return clouds
    
def sampleDisc(clouds, discThick, sbRad = 0.):
    """
    Generate z-axis co-ordinates for particles according to a disc model
    Assumes the galaxy is only weakly triaxial - A,B >> C, but relaxes thin-disc approximation
    
    For strong triaxiality, we would need to include scale height to calculate [x,y,z] from the spaxel [x,y]
    
    ============
    Inputs:
    
    clouds: np.ndarray [x',y',0.,I] output from the rotate clouds routine - which has made thin disc approximation
    
    discThick: float        A constant scale height for the disc
               np.ndarray   The scale height of the disc at radial positions given by sbRad
    
    sbRad (Optional): np.ndarray The radial positions at which discThick is evaluated.
    
    ============
    Outputs:
    
    height: np.ndarray of [z'] co-ordinates
    """
    if isinstance(discThick, (list, tuple, np.ndarray)) or sbRad:
        assert discThick.shape == sbRad.shape, "Disc thickness values and radii dimensions do not match"
    
    rad = np.sqrt(clouds[:,0] ** 2 + clouds[:,1] ** 2)   #radius of each cloud from centre
    if isinstance(discThick, (list, tuple, np.ndarray)):
        interpfunc2 = interpolate.interp1d(sbRad,discThick,kind='linear')
        discThick_here = interpfunc2(rad)
    else:
        discThick_here = discThick
    
    clouds[:,2] = discThick_here * np.random.uniform(-1,1,clouds.shape[0])
    return clouds


def test_skySampler():
    import matplotlib.pyplot as plt
    inc = 75.
    posAng = 240.
    cent = np.array([0.,25.])
    
    cellSize = 0.5
    
    discThick = 0.
    
    rad = np.arange(0.,500.)
    
    nSamps = int(1e6)
    
    # xs, ys, vs, dx,dv,beam,i,
    t0 = time.time()
    model = KinMS.KinMS(100.,100.,100.,cellSize,5.,[1.,1.,0],inc, sbProf=np.where(rad<50.,1.,0.),sbRad=rad,velRad=rad,velProf=rad * 0.,diskThick=discThick,cleanOut=True,
    nSamps=nSamps,posAng=posAng,phaseCen=cent,fileName = 'TestModel.fits')
    t1 = time.time()
    print('Making initial model: ', t1-t0)
    clouds = sampleClouds(model, cellSize, nSamps=nSamps)
    t2 = time.time()
    print('Sampling clouds: ', t2-t1)
    clouds = transformClouds(clouds, cent = cent, inc = inc, posAng = posAng)
    t3 = time.time()
    print('De-rotating clouds ',t3-t2)
    #clouds = sampleDisc(clouds, discThick)
    #print(clouds.shape)
    newMod = KinMS.KinMS(100.,100.,100.,cellSize,5.,[2.,2.,0],inc, velRad=rad,velProf=rad*0.,diskThick=discThick,cleanOut=True,posAng=posAng,phaseCen=cent,
    inClouds = clouds[:,0:3], flux_clouds = clouds[:,3], fileName = 'SSModel.fits')
    t4 = time.time()
    print('Making final model ',t4-t3)
    
    
    fig = plt.figure()
    ax1 = fig.add_subplot(121,aspect='equal')
    ax2 = fig.add_subplot(122,aspect='equal')
    
    Xax = np.arange(-cellSize * model.shape[0]/2, cellSize * model.shape[0]/2, cellSize)
    Yax = np.arange(-cellSize * model.shape[1]/2, cellSize * model.shape[1]/2, cellSize)
    
    levs = np.arange(0,0.0005,0.0001)
    ax1.contourf(Xax,Yax,model.sum(axis=2),cmap="YlOrBr",levels=levs)
    ax1.set_title('Original')
    ax2.contourf(Xax,Yax,newMod.sum(axis=2),cmap="YlOrBr",levels=levs)
    ax2.set_title('SkySampler')
    
    plt.show()

#test_skySampler()
