#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 11:36:51 2024

@author: tom
"""

import electroacPy as ep
import generalToolbox as gtb
import numpy as np

#%% Global parameters
# files
uff_file = "../measurements/zepp_mf_in_box.uff"
step = "../geo/step/zeppelin_mid.step"

# decimate user frequency axis to the one present in uff_file
freq = gtb.freqop.freq_log10(20, 10000, 250)
freq_uff = gtb.freqop.freq_uff(uff_file, freq) 


#%% Mesh CAD data
cad = gtb.meshCAD(step, maxSize=343/10e3/6, minSize=343/10e3/60)
cad.addSurfaceGroup("membrane", [5, 6, 7, 8], 1)
cad.mesh("../geo/mesh/membrane", excludeRemaining=True)

#%% prepare study
infb = ep.loudspeakerSystem(freq_uff)

# set measured points
infb.vibrometry_data("midrange", uff_file, ref2bem=1)
infb.laser_acc["midrange"].plot_pointCloud() # plot current representation without rotation

# change orientation to match CAD mesh
infb.vibrometry_data("midrange", uff_file, rotation=[80, 0, 0], ref2bem=1)
infb.laser_acc["midrange"].plot_pointCloud() # plot updated data


#%% Add observations
infb.study_acousticBEM("infinite_baffle", "../geo/mesh/membrane.msh",
                       "midrange")

infb.observation_polarRadiation("infinite_baffle", "hor_dir", 
                                -90, 270, 5, "y", "x", radius=1,
                                offset=[0, 54.01e-3, 0])
infb.observation_polarRadiation("infinite_baffle", "ver_dir", 
                                -90, 270, 5, "y", "z", radius=1,
                                offset=[0, 54.01e-3, 0])

infb.plot_system("infinite_baffle")

#%% run
infb.run()
ep.save("infinite_baffle", infb)

#%% post-process - horizontal directivity
theta_hor = infb.observation["infinite_baffle"].setup["hor_dir"].theta
pMic_front = infb.get_pMic("infinite_baffle", "hor_dir")
pMic_back = infb.get_pMic("infinite_baffle", "hor_dir")
pMic_hor = pMic_front[:, :73//2+1] + np.flip(pMic_back[:, 73//2:], 1)

gtb.plot.directivityViewer(theta_hor[:73//2+1], freq_uff, pMic_hor)


#%% post-process - vertical directivity
theta_ver = infb.observation["infinite_baffle"].setup["ver_dir"].theta
pMic_front = infb.get_pMic("infinite_baffle", "ver_dir")[:, :73//2+1]
pMic_back = infb.get_pMic("infinite_baffle", "ver_dir")[:, 73//2:]
pMic_ver = pMic_front[:, :73//2+1] + np.flip(pMic_back[:, 73//2::-1], 1)

gtb.plot.directivityViewer(theta_ver[:73//2+1], freq_uff, pMic_ver)


#%% plot at 45 degrees
idx45, _ = gtb.findInArray(theta_hor[:73//2+1], 45)
idx0, _ = gtb.findInArray(theta_hor[:73//2+1], 0)

gtb.plot.FRF(freq_uff, (pMic_hor[:, idx45], 
                        pMic_ver[:, idx45],
                        pMic_hor[:, idx0]), legend=("45 deg - horizontal",
                                                     "45 deg - vertical",
                                                     "on-axis"))







