# KinMS-skySampler
A plugin for Tim Davis' Kinematic Molecular Simulator.
This plugin takes a specified distribution - such as an interferometric data cube, and generates a set of samples that can be input to KinMS to model the data. 
The samples can be drawn according to a uniform, intensity-weighted or custom scheme, and are converted to the plan of the galaxy using a separate routine, so that the sampling algorithm needs no knowledge of the inclination and position angle of the galaxy. This reduces execution time in e.g. an MCMC process.
