#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 15:41:36 2024

This example shows how to use the circuitSolver module. It is based on 
Modified Nodal Analysis (MNA). We simulate both LF as independant acoustic 
radiators playing in the same volume.

The circuitSolver module helps dealing with difficult circuits that are not
part of the predefined setups of the loudspeakerSystem.lem_enclosure() method.
Elements are divided into:
1. components, which generally uses 2 ports,
2. blocks, which are a concatenation of multiple components that form a 
   sub-circuit.

Connections between elements can be either numerical or alphabetical. 

@author: tom
"""

import electroacPy as ep
import generalToolbox as gtb
from math import sqrt, pi

# everything needed for MNA
from electroacPy.circuitSolver.components import acoustic as compa
from electroacPy.circuitSolver.components import electric as compe
from electroacPy.circuitSolver.blocks.electrodynamic import EAD
from electroacPy.circuitSolver.solver import circuit


#%% Initialize system
freq = gtb.freqop.freq_log10(10, 10000, 75)
fst = ep.loudspeakerSystem(freq) # initialize loudspeaker system

#%% Driver definition
# LF/MF
Re = 6
Le = 0.2e-3
Bl = 7.8
Mms = 17.9e-3
Cms = 1.35e-3
Rms = 0.9
Sd = 158e-4
U = sqrt(1.0 * Re) # to get 1W on the woofer

Vs = compe.voltageSource(1, 0, U)
D1 = EAD(1, 0, 2, 3, Le, Re, Cms, Mms, Rms, Bl, Sd)  # LF 1
D2 = EAD(1, 0, 5, 3, Le, Re, Cms, Mms, Rms, Bl, Sd)  # LF 2
R1 = compa.radiator(2, 0, Sd)               # radiation impedance LF 1
R2 = compa.radiator(5, 0, Sd)               # radiation impedance LF 2
B1 = compa.cavity(3, 0, 80e-3)              # box
P1 = compa.port(3, 4, Lp=10e-2, rp=3.5e-2)  # port
RP = compa.radiator(4, 0, pi*3.5e-2**2)     # radiation impedance port


#%% Add LF and port to network 
network_LF = circuit(freq)                      # create network
network_LF.addBlock(D1, D2)                     # add circuit blocks (LF1+LF2)
network_LF.addComponent(Vs, R1, R2, B1, P1, RP) # add components
network_LF.run()

#%% Driver definition
# HF
Re_hf = 3
Le_hf = 0.02e-3
Bl_hf = 3.1
Mms_hf = 0.44e-3
Cms_hf = 117.8e-6
Rms_hf = 0.92
Sd_hf = 9.6e-4

D_HF = EAD(1, 0, 2, 3, Le_hf, Re_hf, Cms_hf, Mms_hf, Rms_hf, Bl_hf, Sd_hf)
R_HF = compa.radiator(2, 0, Sd_hf)
B_HF = compa.cavity(3, 0, 0.02e-3)

#%% Add HF to circuit and run
network_HF = circuit(freq)              # create network
network_HF.addBlock(D_HF)
network_HF.addComponent(Vs, R_HF, B_HF) # because component aren't network dependent, you can re-use already defined component (e.g. Vs)
network_HF.run()

#%% Extract LEM data
# impedance at Vs
Zin_LF = network_LF.getPotential(1)/network_LF.getFlow(1) 
Zin_HF = network_HF.getPotential(1)/network_HF.getFlow(1) 


# LF and port velocities
v_1 = network_LF.getPotential(2) * R1.Gs / Sd      # convert Q to v - LF1
v_2 = network_LF.getPotential(5) * R2.Gs / Sd      # LF2
v_p = network_LF.getPotential(4) * RP.Gs / RP.Sd    # port

# Tweeter velocity
v_hf = network_HF.getPotential(2) * R_HF.Gs / Sd_hf

#%% plot extracted data
gtb.plot.FRF(freq, (Zin_LF, Zin_HF), "abs", 
             legend=("LF - ported", "HF - sealed"),
             title="Input impedance at source U")

gtb.plot.FRF(freq, (v_1*1e3, v_hf*1e3, v_p*1e3), "abs", 
                    legend=("v_lf", "v_hf", "v_p"), 
                    title="driver and port velocities",
                    ylabel="Velocity [mm/s]")


#%% Mesh
# create mesh with separate radiator for each LF
# We use separate LF
cad = gtb.meshCAD("../geo/step/floorstander_piston.step")
cad.addSurfaceGroup("LF1", [12], 1, meshSize=343/5000/6)
cad.addSurfaceGroup("LF2", [14], 2, meshSize=343/5000/6)
cad.addSurfaceGroup("port", [6], 3)
cad.addSurfaceGroup("HF", [13], 4, meshSize=343/5000/6)
cad.addSurfaceGroup("baffle", [9, 10, 11], 5, meshSize=343/2500/6)
cad.mesh("../geo/mesh/floorstander_separate_LF")


#%% Setup simulations
# add "radiator" objects
fst.lem_velocity("LF1", 1, v=v_1)   # add LF velocity to corresponding surface
fst.lem_velocity("LF2", 2, v=v_2)   # add LF velocity to corresponding surface
fst.lem_velocity("port", 3, v=v_p)  # add port velocity
fst.lem_velocity("HF", 4, v=v_hf)   # add tweeter velocity

# setup study 
fst.study_acousticBEM("free-field", 
                      "../geo/mesh/floorstander_separate_LF.msh", 
                      acoustic_radiator =["LF1", "LF2", "port", "HF"])

# setup observations
fst.observation_polarRadiation("free-field", "hor-dir", 
                                -180, 180, 5, "x", "+y", radius=2,
                                offset=[0.2, 0, 0.911])
fst.observation_polarRadiation("free-field", "ver-dir", 
                                -180, 180, 5, "x", "+z", radius=2,
                                offset=[0.2, 0, 0.911])
fst.plot_system("free-field")

#%% Sim
fst.run()
# ep.save("fst_MNA", fst)

#%% Plot
fst.plot_results()


#%% export results - horizontal directivity
import generalToolbox as gtb

pmic_LF1 = fst.get_pMic("free-field", "hor-dir", 1)
pmic_LF2 = fst.get_pMic("free-field", "hor-dir", 2)
pmic_port = fst.get_pMic("free-field", "hor-dir", 3)
pmic_HF = fst.get_pMic("free-field", "hor-dir", 4)
theta = fst.observation["free-field"].setup["hor-dir"].theta

gtb.acoustics.export_directivity("LF1_hor", freq, theta, pmic_LF1)
gtb.acoustics.export_directivity("LF2_hor", freq, theta, pmic_LF2)
gtb.acoustics.export_directivity("HF_hor", freq, theta, pmic_HF)
gtb.acoustics.export_directivity("port_hor", freq, theta, pmic_port)


#%% export results - vertical directivity
pmic_LF1 = fst.get_pMic("free-field", "ver-dir", 1)
pmic_LF2 = fst.get_pMic("free-field", "ver-dir", 2)
pmic_port = fst.get_pMic("free-field", "ver-dir", 3)
pmic_HF = fst.get_pMic("free-field", "ver-dir", 4)
theta = fst.observation["free-field"].setup["ver-dir"].theta


gtb.acoustics.export_directivity("LF1_ver", freq, theta, pmic_LF1)
gtb.acoustics.export_directivity("LF2_ver", freq, theta, pmic_LF2)
gtb.acoustics.export_directivity("HF_ver", freq, theta, pmic_HF)
gtb.acoustics.export_directivity("port_ver", freq, theta, pmic_port)










