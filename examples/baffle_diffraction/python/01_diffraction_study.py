#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 15 10:16:19 2024

@author: tom
"""

import electroacPy as ep
import generalToolbox as gtb


#%% Mesh
fmax = 5e3
Lmax = 343/fmax/6
Lmin = Lmax/10

# cube
cad_cube = gtb.meshCAD("../geo/step/cube.step", Lmin, Lmax)
cad_cube.addSurfaceGroup("piston", [7], 1)
cad_cube.mesh("../geo/msh/cube")

# sphere
cad_sphere = gtb.meshCAD("../geo/step/sphere.step", Lmin, Lmax)
cad_sphere.addSurfaceGroup("piston", [1], 1)
cad_sphere.mesh("../geo/msh/sphere")

# piston
cad_piston = gtb.meshCAD("../geo/step/sphere.step", Lmin, Lmax)
cad_piston.addSurfaceGroup("piston", [1], 1)
cad_piston.mesh("../geo/msh/piston", excludeRemaining=True) # remove everything but surface group 1

#%% setup study
freq = gtb.freqop.freq_log10(20, 5000, 125)
sim = ep.loudspeakerSystem(freq)
sim.lem_velocity("piston", ref2bem=1)

sim.study_acousticBEM("cube", "../geo/msh/cube.msh", "piston")
sim.study_acousticBEM("sphere", "../geo/msh/sphere.msh", "piston")
sim.study_acousticBEM("piston", "../geo/msh/piston.msh", "piston")

sim.observation_polarRadiation(["cube", "sphere", "piston"],
                               "hor-dir", -180, 180, 5, "x", "+y", 
                               offset=[0.1, 0, 0.1])

sim.observation_pressureField(["cube", "sphere", "piston"],
                               "hor-plane", 3, 3, 343/5000/3, "xy",
                               offset=[-1.5, -1.5, 0.1])

# sim.plot_system("cube")
# sim.plot_system("sphere")
# sim.plot_system("piston")

#%% run study
sim.run()

#%% save
ep.save("diffraction_study", sim)

