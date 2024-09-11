"""
ElectroacPy - 2024-04-24_DEV branch

18/07/2024

This script shows how to setup a Lumped-Element network for a bass-reflex enclosure.
The circuit solver class helps to create acoustic circuits that are not pre-defined by the lem_enclosure objects.
"""
import electroacPy as ep
from electroacPy.circuitSolver.circuit import circuit
from electroacPy.circuitSolver.components.acoustic import (cavity, port, open_line, open_line_T,
                                                           closed_line, EAD, radiator)
import generalToolbox as gtb
from generalToolbox.freqop import freq_log10 as f10


mesh = "../geo/mesh_files/606_port.msh"

# DRIVER PARAMETERS
Le, Re, Mms, Cms, Rms, Bl, Sd = 0.16e-3, 3.82, 18.27e-3, 740e-6, 0.47, 6.72, 136.9e-4

# DEFINE COMPONENT IN CIRCUIT
LF = EAD(np0=1,
         np1=2,
         nm=3,
         U=1.0,
         Le=Le,
         Re=Re,
         Cms=Cms,
         Mms=Mms,
         Rms=Rms,
         Bl=Bl,
         Sd=Sd)

RAD1 = radiator(np=2, nm=0, Sd=Sd)         # speaker radiation impedance -> will give Qs, vs
BOX  = cavity(3, 0, Vb=13.53e-3)         # enclosure volume
# PORT = port(3, 4, Lp=18.56e-2, rp=2.5e-2)       # port
PORT = open_line_T(3, 4, 5, 0, Lp=18.56e-2, Sp=3.1415 * 2.5e-2**2)
RAD2 = radiator(5, 0, Sd=3.1415 * 2.5e-2**2)    # radiation from port -> will give Qp, vp

# ASSEMBLE CIRCUIT
freq = f10(10, 2500, 250)
ported_box = circuit(freq)
ported_box.addComponent([LF, RAD1, BOX, PORT, RAD2])

# SOLVE CIRCUIT
ported_box.run()

# GET VELOCITIES TO PUT IN BEM
vs = ported_box.getPotential(2) * RAD1.G / RAD1.Sd   # potential at node 2 (RAD1) * conductance of RAD1 / surface
vp = ported_box.getPotential(5) * RAD2.G / RAD2.Sd   # potential at node 5 (RAD2) * conductance of RAD2 / surface


# BEM STUDY
system = ep.loudspeakerSystem(freq)
system.lem_velocity("LF", ref2bem=1, v=vs)
system.lem_velocity("port", ref2bem=2, v=vp)  # no need to set -vp bc of how the circuit is set up

system.study_exteriorBEM("free-field", mesh, ["LF", "port"])
system.observation_polarRadiation("free-field", "polar_hor",
                                  -180, 180, 5, 'x', '+y',
                                  radius=1.8, offset=[0, 0, 0.175])

system.run()
system.plot_results()
ep.save("saved_data/10_port_circuit_solver", system)


Pin = ported_box.getPotential(3)
P4 = ported_box.getPotential(4)
Qin = (Pin - P4) * PORT.GZt

gtb.acoustics.plot_FRF(freq, Pin/Qin, "dB",
                       ylabel="Input impedance of port (dB)")
