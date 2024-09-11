"""
ElectroacPy - 2024-04-24_DEV branch

18/07/2024

This script shows how set up a vented enclosure LEM/BEM simulation using pre-defined tools.
"""

import electroacPy as ep
import generalToolbox as gtb

# GLOBAL PARAMETERS / FILES
lpm = "../lpm/lpm_lf.txt"    # you can either load Klippel files or set manually each parameter
mesh = "../geo/mesh_files/606_port.msh"
frequency = gtb.freqop.freq_log10(10, 2500, 50)     # start / stop / number of points

# INITIALISATION OF LOUDSPEAKER SYSTEM
system = ep.loudspeakerSystem(frequency)

# LUMPED-ELEMENT PARAMETERS
system.lem_driverFromFile("LF", lpm_data=lpm)
system.lem_enclosure_2("LF_ported",   # name used as a reference for simulations
                       Vb=13.56e-3,         # enclosure volume in m3
                       Lp=18.63e-2,         # port length in m
                       rp=2.5e-2,           # port radius in m
                       setDriver="LF",      # need to set lem_driver into enclosure
                       ref2bem=[1, 2])      # reference to BEM groups, see mesh_from_python file, -> port comes last

# STUDY SETUP
system.study_exteriorBEM("free-field",              # name of study, multiple studies can be defined
                         meshPath=mesh,                   # self-explanatory
                         acoustic_radiator="LF_ported",   # radiators used in simulation, must be referenced to BEM boundaries
                         tol=1e-5)                        # tolerance of GMRES solver, generally doesn't need to be changed

# EVALUATIONS
# polar plot
system.observation_polarRadiation(reference_study="free-field",  # reference to study
                                  observation_name="polar_hor",  # observation name
                                  min_angle=-180,                # start angle
                                  max_angle=180,                 # stop angle
                                  step=5,                        # step
                                  radius=1.8,
                                  on_axis="x",                   # 0 degree facing specified axis
                                  direction="+y",                # direction of rotation, here, +90 degrees will be on "+y" axis
                                  offset=[0, 0, 0.170])          # by default, center of circle at (0, 0, 0)

# field in front of driver
system.observation_pressureField(reference_study="free-field",  # //
                                 observation_name="field_ver",  # //
                                 L1=1,                          # Rectangle length 1
                                 L2=2,                          # Rectangle length 2
                                 step=343/2.5e3/6,              # microphone step (mesh size of rectangular field)
                                 plane="xz",                    # plane on which field is defined
                                 offset=[0.150+0.05, 0, 0.340/2-1])    # by default, rectanle corner at (0, 0, 0)

# field at the back of port
system.observation_pressureField(reference_study="free-field",        # //
                                 observation_name="field_ver_back",   # //
                                 L1=1,                                # //
                                 L2=2,                                # //
                                 step=343/2.5e3/6,                    # //
                                 plane="xz",                          # //
                                 offset=[-1-0.150-0.05, 0, 0.340/2-1])  # //

# display the current acoustic study // might not work or crash on MacOs
system.plot_system("free-field")

# run study
system.run()

# plot results
system.plot_results()

# save results
ep.save("saved_data/01_port", system)

