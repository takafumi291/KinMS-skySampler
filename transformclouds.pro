; SkySampler rotation routines for IDL
FUNCTION transformClouds,inclouds, posAng, inc, cent
 
    ;Calculate the galaxy co-ordinates of clouds from the sky plane. This MUST be used if any of the following conditions are true:
    ;inc != 0
    ;posAng != 90
    ;cent != [0,0]
    ;
    ;This exists as a stand-alone routine since an MCMC fit to the galaxy will likely need to run this every step, 
    ;and the sampleClouds routine is computationally expensive
    ;============
    ;Inputs:
    ;
    ;clouds: np.ndarray The output of the sampleClouds array [x,y,I]
    ;
    ;posAng: 0=<float<360 The position angle of the galaxy major axis, measured from the y-axis of the cube
    ;
    ;inc: 0=<float<90 The inclination of the galaxy, with 90 being edge-on
    ;
    ;cent: [float,float] The photometric centre of the galaxy relative to the centre of the cube
    ;
    ;============
    ;Outputs:
    ;
    ;clouds: np.ndarray Positions and intensities of particles [x',y',I]
    
    clouds =inclouds    ;Make copy of cloud list to avoid modifying original list
    clouds[*,0:1] =[clouds[*,0] - cent[0],clouds[*,1] - cent[1]]    ; Shift clouds from cube position to galaxy-centre position
    posAng = (90-posAng) * !DTOR                                    ; Calculate rotation angle
    xNew = cos(posAng) * clouds[*,0] - sin(posAng) * clouds[*,1]    ; Rotate x-coordinate to align major axis
    yNew = sin(posAng) * clouds[*,0] + cos(posAng) * clouds[*,1]    ; Rotate y-coordinate to align minor axis
    yNew = yNew / cos(inc * !DTOR)                                  ; De-project minor axis
    
    clouds[*,0] = xNew                                              ; Format into an array
    clouds[*,1] = yNew
    return, clouds                                               ; Return the answer
 end
