import electroacPy as ep
import numpy as np
import generalToolbox as gtb


floorstander = ep.load("floorstander_piston")

#%% 
pmic_LF = floorstander.get_pMic("free-field", "hor-dir", [1, 2])
pmic_HF = floorstander.get_pMic("free-field", "hor-dir", 3)
theta = floorstander.observation["free-field"].setup['hor-dir'].theta

gtb.plot.directivityViewer(theta, floorstander.frequency, pmic_LF)


gtb.acoustics.export_directivity("export_LF", floorstander.frequency,
                                 theta, pmic_LF)
gtb.acoustics.export_directivity("export_HF", floorstander.frequency,
                                 theta, pmic_HF)

