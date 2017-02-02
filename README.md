# KinMS-skySampler
A plugin for Tim Davis' Kinematic Molecular Simulator.
This plugin takes a specified sky distribution and converts it into a well-sampled set of cloud positions in the galactic plane, which can be used for KinMS' inClouds mechanism. The input can be a list of positions and intensities (x,y,I), a planar projection of a cube (an array of size XxY containing the intensities along each line of sight) or a cube itself.
