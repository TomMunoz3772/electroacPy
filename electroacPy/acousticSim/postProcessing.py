#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 08:29:39 2024

@author: tom
"""

import numpy as np
from electroacPy.acousticSim.evaluations import evaluations
from electroacPy.acousticSim.bem import bem
from copy import copy, deepcopy
import bempp.api


class postProcess:
    def __init__(self):
        self.TF = {}        # transfer-functions

    def addTransferFunction(self,  name, H, radiatingElement):
        """
        Add a transfer function to the TF dictionnary

        Parameters
        ----------
        name: str
            reference name of the transfer function
        H : ndarray
            transfer function (complex). Must have the same dimension as the 
            BEM/evaluation frequency axis (i.e. same frequency bins)
            
        radiatingElement: int or list of int
            radiating element of the evaluation on which to add the 
            corresponding transfer function.
        Returns
        -------
        None.

        """
        self.TF[name] = {}
        self.TF[name]["H"] = H
        self.TF[name]["radiatingElement"] = radiatingElement
        
        
        
        
        
        
        
        