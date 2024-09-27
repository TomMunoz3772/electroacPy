#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 31 18:58:32 2024

Free-field simulation of a floorstander speaker.
Radiating surfaces are represented as piston.

This is a LF-focused study of example 01_Floorstander_piston. Here we use a
sparser mesh to simulate between 10 to 500 Hz, this allows to speed up a lot 
computations, however, the precision of results at HF will decrease a lot.

In this study, we do not consider the tweeter radiation.

@author: tom
"""

import electroacPy as ep
import generalToolbox as gtb
from math import sqrt

#%% Initialize system
freq = gtb.freqop.freq_log10(10, 10000, 50)
floorstander = ep.loudspeakerSystem(freq)

#%% Import CAD and create mesh
max_s = 343/500/6 
min_s = max_s/10

step_file = "../geo/step/floorstander_piston.step" 
cad = gtb.meshCAD(step_file, min_s, max_s)
cad.addSurfaceGroup("LF", [12, 14], 1)
cad.addSurfaceGroup("port", [6], 2)
cad.addSurfaceGroup("HF", [13], 3)
cad.mesh("../geo/mesh/floorstander_piston_LF")


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
floorstander.lem_driver("MW19P", U, Le, Re, Cms, Mms, Rms, Bl, Sd)
floorstander.lem_enclosure("ported_box", 80e-3, Lp=10e-2, rp=3.5e-2, Nd=2, 
                    setDriver="MW19P", ref2bem=[1, 2])

#%% Observation setup 
floorstander.study_acousticBEM("free-field", 
                               "../geo/mesh/floorstander_piston_LF.msh", 
                               "ported_box")

floorstander.observation_polarRadiation("free-field", "hor-dir", 
                               -180, 180, 5, "x", "+y", radius=2,
                               offset=[0, 0, 0.91])
floorstander.observation_polarRadiation("free-field", "ver-dir", 
                               -180, 180, 5, "x", "+z", radius=2,
                               offset=[0, 0, 0.91])

floorstander.observation_pressureField("free-field", "hor-pf", 10, 10, 
                                       343/250/6, "xy",
                                       offset=[-5, -5, 0.91])
floorstander.observation_pressureField("free-field", "ver-pf", 10, 10, 
                                       343/250/6, "xz",
                                       offset=[-5, 0, -5+0.6])

floorstander.plot_system("free-field")

#%% Sim
floorstander.run()
ep.save("floorstander_piston_LF", floorstander)

#%% Plot
pmic = floorstander.get_pMic("free-field", "hor-dir", [1, 2])

gtb.plot.FRF(freq, (pmic[:, 73//2], pmic[:, 0],
                    pmic[:, 50]), xlim=(10, 500),
             legend=("on-axis", "-180 degrees", "+50 degrees"),
             title='Range of study where mesh is "correct"')

gtb.plot.FRF(freq, (pmic[:, 73//2], pmic[:, 0],
                    pmic[:, 50]), xlim=(500, 10000), 
             legend=("on-axis", "-180 degrees", "+50 degrees"),
             title='Range of study where mesh is too large')

#%% Directivity
floorstander.plot_results(observation="hor-dir")


