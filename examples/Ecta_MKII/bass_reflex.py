#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 12:12:22 2022

@author: tom
"""

import electroacPy as ep
from generalToolbox.freqop import freq_log10
from numpy import pi


freq = freq_log10(20, 2500, 75)
ecta = ep.loudspeakerSystem(freq)  # freq is optional, default being freq_log10(20, 2500, 50)

## midrange
Re = 5.9
Le = 0.13e-3
Bl = 6.1
Mms = 4.7e-3
Cms = 1.24e-3
Rms = 0.54
Sd = 58e-4
ecta.lem_driver("midrange", 1, Le, Re, Cms, Mms, Rms, Bl, Sd)

# before setting enclosure volume:
# ecta.driver['midrange'].sealedAlignment()
Vb_mr = 4.9e-3  # volume box mid-range // Qtc = 0.43
ecta.lem_enclosure("mr_sealed", Vb_mr, setDriver="midrange", ref2bem=4)

## woofer
Re = 3.3
Le = 0.5e-3
Bl = 5.9
Mms = 19.1e-3
Cms = 1.2e-3
Rms = 1.1
Sd = 154e-4
ecta.lem_driver("woofer", 1, Le, Re, Cms, Mms, Rms, Bl, Sd)

# ported alignment estimator
# ecta.driver['woofer'].portedAlignment()
Vb_woofer = 30e-3
Lp = 160e-3
Sp = pi*(68e-3 / 2)**2

ecta.lem_enclosure("wf_ported", Vb_woofer, Lp, Sp, setDriver="woofer", ref2bem=[3, 6])

## acoustic sim
msh = "ecta_mesh_500Hz.msh"  # mesh can be opened as a regular .txt file - useful for ref2bem


ecta.study_exteriorBEM('free-field', msh,
                       ["mr_sealed", "wf_ported"])
ecta.observation_polarRadiation('free-field', 'directivity',
                                -180, 180, 1, 'xy',
                                radius=1.8, offset=[0.200, 0.110, 1])
ecta.observation_pressureField('free-field', 'v_propag',
                               4, 3, 343/1000/5, 'xz', offset=[-1, 0.110, -1])


# ecta.plot_system('free-field')
ecta.run()

## save and plot results
ep.save("ecta_500Hz", ecta)
ecta.plot_results()
# depending on computer OS, PyVista might crash when closing pressureField evaluations
# (not encountered on Windows, but seems to happen a lot on MacOs) -> best to save before plotting results


## Let's add crossovers
ecta.filter_network('LF', ref2bem=[3, 6], ref2study='free-field') # applied to woofer and port
ecta.filter_addLowPass('LF', 'lp1', 1, fc=2*ecta.driver['midrange'].Fs)

ecta.filter_network('MF', ref2bem=4, ref2study='free-field')
ecta.filter_addHighPass('MF', 'hp1', 1, fc=2*ecta.driver['midrange'].Fs)
ecta.filter_addGain('MF', 'db1', -5)
ecta.filter_addPhaseFlip('MF', 'pi')


## Some new plots
ecta.plot_xovers(['LF', 'MF'])
ecta.plot_results(observation='directivity')
ecta.plot_results(observation='directivity', bypass_xover=True, radiatingSurface=4)


## Load results
ecta_II = ep.load("ecta_500Hz")  # ecta_II will behave as any regular loudspeakerSystem() object
ecta_II.export_results("ecta_500Hz - CSV EXPORT")
