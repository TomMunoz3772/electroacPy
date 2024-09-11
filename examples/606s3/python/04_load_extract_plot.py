"""
ElectroacPy - 2024-04-24_DEV branch

18/07/2024

This script shows how load saved simulation, extract and plot microphone pressure
"""

import electroacPy as ep
import generalToolbox as gtb


# LOAD DATA
system_ported = ep.load("saved_data/01_port")
system_abr    = ep.load("saved_data/02_ABR")
system_sealed = ep.load("saved_data/03_sealed")
frequency = system_sealed.frequency

# EXTRACT TOTAL MICROPHONE PRESSURE OF THE POLAR RADIATION
pMic_ported = system_ported.get_pMic("free-field", "polar_hor")
pMic_abr    = system_abr.get_pMic("free-field", "polar_hor")
pMic_sealed = system_sealed.get_pMic("free-field", "polar_hor")

# PLOT USING GENERAL TOOLBOX
gtb.acoustics.plot_FRF(frequency, H=(pMic_ported[:, 73//2],  # [:, 73//2] get SPL for all frequencies at index 36 (theta=0 in studies)
                                     pMic_abr[:, 73//2],
                                     pMic_sealed[:, 73//2]),
                       transformation="SPL",
                       labels=("ported", "ABR", "sealed"))

# EXTRACT AND PLOT SPECIFIC SURFACES
PORT = system_ported.get_pMic("free-field", "polar_hor", radiatingSurface=2)
ABR  = system_abr.get_pMic("free-field", "polar_hor", radiatingSurface=2)

gtb.acoustics.plot_FRF(frequency, (PORT[:, 73//2], ABR[:, 73//2]),
                       labels=("port pressure", "ABR pressure"))


