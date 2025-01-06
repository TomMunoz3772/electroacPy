#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 08:29:39 2024

@author: tom
"""

class postProcess:
    def __init__(self):
        """
        Store transfer function of corresponding radiating elements.
        Should be set in a evaluation.plot() method call.
        
        ex:
        EV = evaluationObject
        EV.plot(evaluations=[], radiatingElement=[1, 2, 3], 
                processing=postProcessObj)
        
        -> will plot all evaluations with summed elements' TF for 1, 2 and 3.
        """
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
        
        
        
        
        
        
        
        