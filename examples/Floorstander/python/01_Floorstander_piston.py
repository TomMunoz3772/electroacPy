#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 31 18:58:32 2024

Free-field simulation of a floorstander speaker.
Radiating surfaces are represented as piston.

@author: tom
"""

import electroacPy as ep
import generalToolbox as gtb
from math import sqrt

#%% Initialize system
freq = gtb.freqop.freq_log10(10, 10000, 50)
floorstander = ep.loudspeakerSystem(freq)

#%% Import CAD and create mesh
max_s = 343/1000/6 
min_s = max_s/10

step_file = "../geo/step/floorstander_piston.step" 
cad = gtb.meshCAD(step_file, min_s, max_s)
cad.addSurfaceGroup("LF", [12, 14], 1, meshSize=343/5000/6)
cad.addSurfaceGroup("port", [6], 2)
cad.addSurfaceGroup("HF", [13], 3, meshSize=343/5000/6)
cad.addSurfaceGroup("baffle", [9, 10, 11],
                    4, meshSize=343/3500/6)
cad.mesh("../geo/mesh/floorstander_piston")


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

# HF
Re = 3
Le = 0.02e-3
Bl = 3.1
Mms = 0.44e-3
Cms = 117.8e-6
Rms = 0.92
Sd = 9.6e-4
floorstander.lem_driver("TW29RN", U, Le, Re, Cms, Mms, Rms, Bl, Sd)
floorstander.lem_enclosure("sealed_tw", 0.0122e-3, setDriver="TW29RN", ref2bem=3)

#%% Observation setup 
floorstander.study_acousticBEM("free-field", "../geo/mesh/floorstander_piston.msh", 
                      ["ported_box", "sealed_tw"])

floorstander.observation_polarRadiation("free-field", "hor-dir", 
                               -180, 180, 5, "x", "+y", radius=2,
                               offset=[0, 0, 0.91])
floorstander.observation_polarRadiation("free-field", "ver-dir", 
                               -180, 180, 5, "x", "+z", radius=2,
                               offset=[0, 0, 0.91])

floorstander.observation_pressureField("free-field", "hor-pf", 10, 10, 343/250/6, "xy",
                              offset=[-5, -5, 0.91])
floorstander.observation_pressureField("free-field", "ver-pf", 10, 10, 343/250/6, "xz",
                              offset=[-5, 0, -5+0.6])

floorstander.plot_system("free-field")

#%% Sim
floorstander.run()
# ep.save("floorstander_piston", floorstander)

#%% Plot
floorstander.plot_results(radiatingElement=[1, 2])

#%% Export data - Horizontal directivity
pmic_LF = floorstander.get_pMic("free-field", "hor-dir", [1, 2])
pmic_HF = floorstander.get_pMic("free-field", "hor-dir", 3)
theta = floorstander.observation["free-field"].setup['hor-dir'].theta

gtb.acoustics.export_directivity("export_LF_hor", floorstander.frequency,
                                 theta, pmic_LF)
gtb.acoustics.export_directivity("export_HF_hor", floorstander.frequency,
                                 theta, pmic_HF)

#%% Export data - Vertical directivity
pmic_LF = floorstander.get_pMic("free-field", "ver-dir", [1, 2])
pmic_HF = floorstander.get_pMic("free-field", "ver-dir", 3)
theta = floorstander.observation["free-field"].setup['ver-dir'].theta

gtb.acoustics.export_directivity("export_LF_ver", floorstander.frequency,
                                 theta, pmic_LF)
gtb.acoustics.export_directivity("export_HF_ver", floorstander.frequency,
                                 theta, pmic_HF)


