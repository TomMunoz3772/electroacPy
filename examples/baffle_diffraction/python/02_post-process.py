#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 15 10:46:55 2024

@author: tom
"""

import electroacPy as ep
import generalToolbox as gtb
import numpy as np

#%% load and extract results from study
study = ep.load("diffraction_study")
cube = study.observation["cube"].setup["hor-dir"].pMic[:, 73//2] 
sphere = study.observation["sphere"].setup["hor-dir"].pMic[:, 73//2] 
piston = study.observation["piston"].setup["hor-dir"].pMic[:, 73//2] 
piston_mirror = study.observation["piston"].setup["hor-dir"].pMic[:, 0]
piston_inf = piston + piston_mirror

#%% plot
gtb.plot.FRF(study.frequency, 
             (cube/piston_inf, sphere/piston_inf, piston_inf/piston_inf), 
             transformation="dB", 
             legend=("cube", "sphere", "infinite baffle"),
             title="estimated baffle diffraction",
             ylim=(-12, 12), 
             xticks=(20, 50, 100, 500, 1000, 5000),
             yticks=np.arange(-12, 15, 3))


