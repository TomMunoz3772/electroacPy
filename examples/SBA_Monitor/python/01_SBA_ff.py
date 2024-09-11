#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 31 18:58:32 2024

Free-field simulation of a bookshelf speaker.
Radiating surfaces are represented as piston.

Info and data: http://www.troelsgravesen.dk/MW19P-8.htm

@author: tom
"""

import electroacPy as ep
import generalToolbox as gtb
from math import sqrt

#%% Initialize system
freq = gtb.freqop.freq_log10(10, 10000, 50)
sba = ep.loudspeakerSystem(freq)

#%% Import CAD and create mesh
max_s = 343/1000/6 
min_s = max_s/10

step_file = "geo/step/SBA_M_ext.step" 
cad = gtb.meshCAD(step_file, min_s, max_s)
cad.addSurfaceGroup("LF", [15], 1, meshSize=343/8000/6)
cad.addSurfaceGroup("port", [16], 2)
cad.addSurfaceGroup("HF", [14], 3, meshSize=343/8000/6)
cad.addSurfaceGroup("baffle", [7, 11, 12, 13, 5, 3, 10, 9],
                    4, meshSize=343/3500/6)
cad.mesh("geo/mesh/SBA_M_ext")


#%% Driver definition
# LF/MF
Re = 6
Le = 0.2e-3
Bl = 7.8
Mms = 17.9e-3
Cms = 1.35e-3
Rms = 0.9
Sd = 158e-4
U = sqrt(1.0 * Re)
sba.lem_driver("MW19P", U, Le, Re, Cms, Mms, Rms, Bl, Sd)
sba.lem_enclosure_2("ported_box", 24e-3, Lp=16e-2, rp=2.5e-2, 
                    setDriver="MW19P", ref2bem=[1, 2])

# HF
Re = 3
Le = 0.02e-3
Bl = 3.1
Mms = 0.44e-3
Cms = 117.8e-6
Rms = 0.92
Sd = 9.6e-4
sba.lem_driver("TW29RN", U, Le, Re, Cms, Mms, Rms, Bl, Sd)
sba.lem_enclosure_2("sealed_tw", 0.0122e-3, setDriver="TW29RN", ref2bem=3)

#%% Observation setup 
sba.study_exteriorBEM("free-field", "geo/mesh/SBA_M_ext.msh", 
                      ["ported_box", "sealed_tw"])

sba.observation_polarRadiation("free-field", "hor-dir", 
                               -180, 180, 5, "x", "+y", radius=2,
                               offset=[280e-3/2, 260e-3/2, 460e-3/2])
sba.observation_polarRadiation("free-field", "ver-dir", 
                               -180, 180, 5, "x", "+z", radius=2,
                               offset=[280e-3/2, 260e-3/2, 460e-3/2])

sba.observation_pressureField("free-field", "hor-pf", 4, 4, 343/750/6, "xy",
                              offset=[-2+140e-3, -2+130e-3, 230e-3])
sba.observation_pressureField("free-field", "ver-pf", 4, 4, 343/750/6, "xz",
                              offset=[-2+140e-3, 130e-3, -2+230e-3])


#%% Sim
sba.run()
# ep.save("SBA_M", sba)

#%% Plot
sba.plot_results(radiatingSurface=[1, 2])
