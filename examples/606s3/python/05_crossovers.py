"""
ElectroacPy - 2024-04-24_DEV branch

18/07/2024

This script simulate the radiation of a 606s3 LF and tweeter, crossovers are done using the xover class.
Because the crossovers are applied through transfer function, the change of impedance at LF due to using passive
component is not taken into account; however, digital filters are true to their response.

Filters are generally added after simulation. This means that the BEM solver can be run only once.
"""
import electroacPy as ep
import generalToolbox as gtb
import numpy as np

# INITIAL MESH
cad = gtb.meshCAD("../geo/step_files/606_PORT_TWEETER.step", maxSize=343 / 4e3 / 6, minSize=343 / 4e3 / 60)
cad.addSurfaceGroup("LF", [9], 1)
cad.addSurfaceGroup("port", [7], 2)
cad.addSurfaceGroup("tweeter", [8], 3)
cad.mesh("../geo/mesh_files/606_port_tweeter")

# GLOBAL PARAMETERS / FILES
lpm_LF = "../lpm/lpm_lf.txt"  # you can either load Klippel files or set manually each parameter
lpm_HF = "../lpm/lpm_hf.txt"
mesh = "../geo/mesh_files/606_port_tweeter.msh"
frequency = gtb.freqop.freq_log10(10, 7500, 75)  # start / stop / number of points

# INITIALISATION OF LOUDSPEAKER SYSTEM
system = ep.loudspeakerSystem(frequency)

# LUMPED-ELEMENT PARAMETERS
system.lem_driverFromFile("LF", lpm_data=lpm_LF)
system.lem_driverFromFile("HF", lpm_data=lpm_HF)
system.lem_enclosure_2("LF_ported",  # name used as a reference for simulations
                       Vb=13.56e-3,  # enclosure volume in m3
                       Lp=18.63e-2,  # port length in m
                       rp=2.5e-2,  # port radius in m
                       setDriver="LF",  # need to set lem_driver into enclosure
                       ref2bem=[1, 2])  # reference to BEM groups, see mesh_from_python file, -> port comes last

system.lem_enclosure_2("HF_sealed", 0.007e-3, setDriver="HF", ref2bem=3)

# STUDY SETUP
system.study_exteriorBEM("free-field",  # name of study, multiple studies can be defined
                         meshPath=mesh,  # self-explanatory
                         acoustic_radiator=["LF_ported", "HF_sealed"],
                         # radiators used in simulation, must be referenced to BEM boundaries
                         tol=1e-5)  # tolerance of GMRES solver, generally doesn't need to be changed

# EVALUATIONS
# polar plot
system.observation_polarRadiation(reference_study="free-field",  # reference to study
                                  observation_name="polar_hor",  # observation name
                                  min_angle=-180,  # start angle
                                  max_angle=180,  # stop angle
                                  step=5,  # step
                                  radius=1.8,
                                  on_axis="x",  # 0 degree facing specified axis
                                  direction="+y",  # direction of rotation, here, +90 degrees will be on "+y" axis
                                  offset=[0, 0, 0.170])  # by default, center of circle at (0, 0, 0)

system.run()

# get listening window
# all drivers
pMic_no_filter = system.get_pMic("free-field", "polar_hor")
pMic_lst_no_filter = np.mean(pMic_no_filter[:, 73 // 2 - 6:73 // 2 + 7], 1)  #  average over +-30 degrees

# tweeter only
pMic_tweeter = system.get_pMic("free-field", "polar_hor", radiatingSurface=3)
pMic_tweeter_lst = np.mean(pMic_tweeter[:, 73//2-6:73//2+7], 1)

# LF + port only
pMic_LF = system.get_pMic("free-field", "polar_hor", radiatingSurface=[1, 2])
pMic_LF_lst = np.mean(pMic_LF[:, 73//2-6:73//2+7], 1)

# plot
gtb.acoustics.plot_FRF(frequency, (pMic_lst_no_filter,
                                   pMic_tweeter_lst,
                                   pMic_LF_lst), labels=("total - no filter", "tweeter", "LF"))


# filter
# from the previous plot, the midrange radiate an average of 4 dB more than the tweeter (mostly due to baffle effect),
# we set the crossover at 800 Hz
system.filter_network("xo_LF", ref2bem=[1, 2], ref2study="free-field")
system.filter_addLowPass("xo_LF", "lp1", 1, 800)

system.filter_network("xo_HF", ref2bem=3, ref2study="free-field")
system.filter_addHighPass("xo_HF", "hp1", 1, 800)
system.filter_addDelay("xo_HF", "dt1", 0.02e-3)

# get listening window
# all drivers
pMic = system.get_pMic("free-field", "polar_hor")
pMic_lst = np.mean(pMic[:, 73 // 2 - 6:73 // 2 + 7], 1)  #  average over +-30 degrees

# plot updated pressure
gtb.acoustics.plot_FRF(frequency, (pMic_lst_no_filter, pMic_lst),
                       labels=("total - no filter", "total - filter"))

