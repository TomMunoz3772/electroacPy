#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 11:00:55 2024

@author: tom
"""

import electroacPy as ep
import generalToolbox as gtb

#%% Load data
piston = ep.load("floorstander_piston")
piston_LF = ep.load("floorstander_piston_LF")
cone = ep.load("floorstander_cone")


#%% extract pressure from the horizontal directivity
pmic_pst = piston.get_pMic("free-field", "hor-dir", radiatingElement=1)         # piston
pmic_pst_LF = piston_LF.get_pMic("free-field", "hor-dir", radiatingElement=1)   # piston LF
pmic_cne = cone.get_pMic("free-field", "hor-dir", radiatingElement=1)           # cone

pmic_pst_p = piston.get_pMic("free-field", "hor-dir", radiatingElement=2)       # piston
pmic_pst_LF_p = piston_LF.get_pMic("free-field", "hor-dir", radiatingElement=2) # piston LF
pmic_cne_p = cone.get_pMic("free-field", "hor-dir", radiatingElement=2)         # cone


#%% Plot
# the cone pressure is normalized down to the actual effective radiating sur-
# face from Thiele/Small parameters.
gtb.plot.FRF(piston.frequency, 
            (pmic_pst[:, 73//2], 
             pmic_cne[:, 73//2],
             pmic_pst_p[:, 73//2],
             pmic_cne_p[:, 73//2]), 
            legend=("piston - driver", "cones -driver",
                    "piston - port", "'cone' - port"),
            ylim=(10, 110), xlim=(10, 10000),
            title="individual radiating elements")

gtb.plot.FRF(piston.frequency,
             (pmic_pst[:, 73//2]+ pmic_pst_p[:, 73//2],
             pmic_cne[:, 73//2] + pmic_cne_p[:, 73//2]),
             legend=("pistons", "cones"), 
             ylim=(10, 110), xlim=(20, 10000), 
             title="driver + port")

gtb.plot.FRF(piston.frequency,
             (pmic_pst[:, 0]+pmic_pst_p[:, 0], 
              pmic_pst_LF[:, 0]+pmic_pst_LF_p[:, 0]),
             legend=("thin mesh", "large mesh"),
             title="Off-axis precision difference between large and thin mesh")
