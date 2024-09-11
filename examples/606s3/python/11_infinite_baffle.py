"""
ElectroacPy - 2024-04-24_DEV branch

18/07/2024

This script shows how to return the infinite baffle pressure. As the specific Green function is not defined
in bempp-cl, we use the back radiation of a flat shell as the imaginary source pressure.

Because we use an open shell (and not a completely closed, water-tight mesh), we can greatly reduce the element
size without increasing too much the computation time.
"""

import electroacPy as ep
import generalToolbox as gtb
import numpy as np

# MESH
cad = gtb.meshCAD("../geo/step_files/speaker_cone.step", minSize=343/8e3/60, maxSize=343/8e3/6)
cad.addSurfaceGroup("membrane", [1, 2, 3, 4, 5], 1)
cad.mesh("../geo/mesh_files/speaker_cone", excludeRemaining=True)


# lpm and mesh
lpm = "../lpm/lpm_lf.txt"
# mesh = "../geo/mesh_files/speaker_cone.msh"
mesh = "../geo/mesh_files/606_baffle.msh"

# BEM
frequency = gtb.freqop.freq_log10(20, 12e3, 125)
system = ep.loudspeakerSystem(frequency)
system.lem_driverFromFile("LF", lpm)
system.lem_enclosure_2("LF_sealed", 13.53e-3, ref2bem=1, setDriver="LF")

system.study_exteriorBEM("inf_baffle", mesh,
                         "LF_sealed")

system.observation_polarRadiation("inf_baffle", "directivity",
                                  -180, 180, 5, "x", "+y", offset=[0.15, 0, 0.17])

system.observation_pressureField("inf_baffle", "hor_plane", 2, 2, 343/2500/6,
                                 'xy', offset=[-1, -1, 0.17])

# system.plot_system("inf_baffle")
system.run()


# because the simulation is done in free-field, we need to sum the pressure behind and in-front of the imaginary
# infinite baffle
pMic = system.get_pMic("inf_baffle", "directivity")
pMic_front = pMic[:, 73//4:3*73//4+1]  # from -90 to +90 -> array starts at -180 with a 5 degrees step
pMic_back = np.concatenate((pMic[:, 73//4::-1], pMic[:, -1:-73//4:-1]), 1)  # we get back pressure from +90 to -90
pMic_inf_baffle = pMic_front+pMic_back

gtb.acoustics.plot_directivity(frequency, np.arange(-90, 95, 5), pMic_inf_baffle.T,
                               ylabel="Angle [degrees]", levels=np.arange(35, 75, 3))

gtb.acoustics.plot_FRF(frequency, pMic_inf_baffle[:, 37//2])
