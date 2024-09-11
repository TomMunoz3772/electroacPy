#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 11:00:55 2024

@author: tom
"""

import electroacPy as ep
import generalToolbox as gtb

#%% Load data
sba_piston = ep.load("SBA_M")
sba_cone = ep.load("SBA_M_cones")


#%% extract pressure from the horizontal directivity
pmic_pst = sba_piston.get_pMic("free-field", "hor-dir", radiatingSurface=1)     # piston
pmic_cne = sba_cone.get_pMic("free-field", "hor-dir", radiatingSurface=1)       # cone

pmic_pst_p = sba_piston.get_pMic("free-field", "hor-dir", radiatingSurface=2)   # piston
pmic_cne_p = sba_cone.get_pMic("free-field", "hor-dir", radiatingSurface=2)     # cone


#%% Plot
# the cone pressure is normalized down to the actual effective radiating sur-
# face from Thiele/Small parameters.
gtb.acoustics.plot_FRF(sba_piston.frequency, 
                       (pmic_pst[:, 73//2], 
                        pmic_cne[:, 73//2] * 10**(-4.459/20),
                        pmic_pst_p[:, 73//2],
                        pmic_cne_p[:, 73//2],
                        pmic_pst[:, 73//2]+pmic_pst_p[:, 73//2],
                        pmic_cne[:, 73//2]* 10**(-4.459/20)
                        +pmic_cne_p[:, 73//2]), 
                       labels=("piston - driver", "cones -driver",
                               "piston - port", "'cone' - port", 
                               "sum - piston", "sum - cone"),
                       ylim=(10, 110), xlim=(20, 10000))
