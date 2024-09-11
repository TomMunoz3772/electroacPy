#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 15:41:36 2024

@author: tom
"""

import electroacPy as ep
import generalToolbox as gtb
from math import sqrt, pi
from electroacPy.circuitSolver.components import acoustic as compa
from electroacPy.circuitSolver.components import electric as compe
from electroacPy.circuitSolver.blocks.electrodynamic import EAD
from electroacPy.circuitSolver.blocks.electric import lowpass_butter, highpass_butter
from electroacPy.circuitSolver.solver import circuit


#%% Initialize system
freq = gtb.freqop.freq_log10(10, 10000, 75)
sba = ep.loudspeakerSystem(freq) # initialize loudspeaker system

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

Vs = compe.voltageSource("A", 0, U)
LP = lowpass_butter("A", "B", 3, 2000, Re)
LF = EAD("B", 0, "D", "E", Le, Re, Cms, Mms, Rms, Bl, Sd)
rad_lf = compa.radiator("D", 0, Sd)
cab_lf = compa.cavity("E", 0, 24e-3)
port_lf = compa.port("E", "F", 16.1e-2, 2.5e-2)
rad_p = compa.radiator("F", 0, Sd=pi*2.5e-2**2)


#%% Add LF and port to network 
network = circuit(freq)                                  # create network
network.addBlock(LP, LF)                                 # add circuit blocks
network.addComponent(Vs, rad_lf, cab_lf, port_lf, rad_p) # add components


#%% Driver definition
# HF
Re_hf = 3
Le_hf = 0.02e-3
Bl_hf = 3.1
Mms_hf = 0.44e-3
Cms_hf = 117.8e-6
Rms_hf = 0.92
Sd_hf = 9.6e-4

HP = highpass_butter("A", "G", 3, 2000, Re_hf)
HF = EAD("G", 0, "I", "J", Le_hf, Re_hf, Cms_hf, Mms_hf, Rms_hf, Bl_hf, Sd_hf)
rad_hf = compa.radiator("I", 0, Sd_hf)
cab_hf = compa.cavity("J", 0, 0.02e-3)

#%% Add HF to circuit and run
network.addBlock(HP, HF)
network.addComponent(rad_hf, cab_hf)
network.run()

#%% Extract LEM data
# impedance at Vs
Zin = network.getPotential("A")/network.getFlow("A")

# LF and port velocities
v_lf = network.getPotential("D") * rad_lf.Gs / Sd      # convert Q to v
v_p = network.getPotential("F") * rad_p.Gs / rad_p.Sd

# Tweeter velocity
v_hf = network.getPotential("I") * rad_hf.Gs / Sd_hf

# crossover's transfer function
H_lp = network.getPotential("B") / network.getPotential("A")
H_hp = network.getPotential("G") / network.getPotential("A")


#%% plot extracted data
gtb.acoustics.plot_FRF(freq, Zin, "abs", 
                       title="Input impedance at source U")
gtb.acoustics.plot_FRF(freq, (H_lp, H_hp, H_lp+H_hp), 
                       legend=("LP", "HP", "sum"), 
                       title="crossover's transfer function")
gtb.acoustics.plot_FRF(freq, (v_lf*1e3, v_hf*1e3, v_p*1e3), "abs", 
                        legend=("v_lf", "v_hf", "v_p"), 
                        title="driver and port velocities",
                        ylabel="Velocity [mm/s]")

#%% Setup simulations
# add "radiator" objects
sba.lem_velocity("LF", 1, v=v_lf)   # add LF velocity to corresponding surface
sba.lem_velocity("port", 2, v=v_p)  # add port velocity
sba.lem_velocity("HF", 3, v=v_hf)   # add tweeter velocity

# setup study 
sba.study_exteriorBEM("free-field", "../geo/mesh/SBA_M_ext_cones.msh", 
                      acoustic_radiator =["LF", "port", "HF"])

# setup observations
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

sba.plot_system("free-field")

#%% Sim
sba.run()
# ep.save("SBA_MNA", sba)

#%% Plot
sba.plot_results()
