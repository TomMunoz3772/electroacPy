#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 15 10:46:55 2024

@author: tom
"""

import electroacPy as ep
import generalToolbox as gtb


#%% load and extract results from study
study = ep.load("diffraction_study")
cube = study.get_pMic("cube", "hor-dir")[:, 73//2]         # on-axis
sphere = study.get_pMic("sphere", "hor-dir")[:, 73//2]     # ""
piston = study.get_pMic("piston", "hor-dir")[:, 73//2]     # ""
piston_mirror = study.get_pMic("piston", "hor-dir")[:, 0]  # off-axis (180°)
piston_tot = piston + piston_mirror                        # equivalent infinite baffle


#%% plot
gtb.acoustics.plot_FRF(study.frequency, (cube/piston_tot, sphere/piston_tot,
                                         piston_tot/piston_tot), 
                       transformation="dB", legend=("cube", "sphere", 
                                                    "infinite baffle"),
                       title="impact of enclosure shape on loudspeaker gain")


