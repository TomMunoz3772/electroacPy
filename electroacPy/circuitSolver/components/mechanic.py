#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 10:40:49 2024

@author: tom
"""

from numpy import array, zeros, ones, repeat
from numpy import max as npmax
from generalToolbox.freqop import laplace
from generalToolbox import parallel

class voltageSource:
    def __init__(self, np, nm, value):
        """
        Create a resistor component.

        Parameters
        ----------
        np : int
            positive node.
        nm : int
            negative node.
        value : float
            resustance value.

        Returns
        -------
        None.

        """
        
        self.np = np
        self.nm = nm
        self.value = value
        self.G = 1
        self.Gs = None
        
        # create stamp and relative informations
        self.stamp_G = array([[0, 0, 1], [0, 0, -1], [1, -1, 0]])
        self.stamp_I = array([[0], [0], [value]])
        self.contribute = ["G", "Is"]  
        self.vsource = 1
    
        
    def init_component(self, frequency):
        self.Gs = self.G * ones(len(frequency)) 
        
    
    def update_stamp(self, node_list, M, nbsource):
        """
        Update the component's stamp in a matrix of the system's size.

        Parameters
        ----------
        node_list : list
            list of nodes.

        Returns
        -------
        None.

        """
        np = self.np
        nm = self.nm
        maxNode = npmax(node_list)
        
        # G stamp
        self.stamp_G = zeros([maxNode+M, maxNode+M], dtype=complex)
        
        if self.np != 0:
            self.stamp_G[self.np-1, maxNode+nbsource] = 1
            self.stamp_G[maxNode+nbsource, self.np-1] = 1
        if self.nm != 0:
            self.stamp_G[maxNode+nbsource, self.nm-1] = -1        
        
        # I stamp
        self.stamp_I = zeros([maxNode+M, 1], dtype=complex) 
        self.stamp_I[maxNode+nbsource] = 1
        

class resistance:
    def __init__(self, np, nm, value):
        """
        Create a resistor component.

        Parameters
        ----------
        np : int
            positive node.
        nm : int
            negative node.
        value : float
            resustance value.

        Returns
        -------
        None.

        """
        
        self.np = np
        self.nm = nm
        self.value = value
        self.G = 1/value
        self.Gs = None
        
        # create stamp and relative informations
        self.stamp_G = array([[1, -1], [-1, 1]])
        self.contribute = ["G"]        
        self.vsource = 0

    def init_component(self, frequency):
        """
        initialize resistor within circuit

        Parameters
        ----------
        frequency : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        self.Gs = self.G * ones(len(frequency))  
        
    def update_stamp(self, node_list, M):
        """
        Update the component's stamp in a matrix of the system's size.
 
        Parameters
        ----------
        node_list : list
            list of nodes.
 
        Returns
        -------
        None.
 
        """
        np = self.np
        nm = self.nm
        maxNode = npmax(node_list)
        self.stamp_G = zeros([maxNode+M, maxNode+M], dtype=complex)
        
        if np != 0:
            self.stamp_G[np - 1, np - 1] = 1
        if nm != 0:
            self.stamp_G[nm - 1, nm - 1] = 1
        if np != 0 and nm != 0:
            self.stamp_G[np - 1, nm - 1] = -1
            self.stamp_G[nm - 1, np - 1] = -1


class inductance:
    def __init__(self, np, nm, value):
        """
        Create a inductor component.

        Parameters
        ----------
        np : int
            positive node.
        nm : int
            negative node.
        value : float
            resustance value.

        Returns
        -------
        None.

        """
        
        self.np = np
        self.nm = nm
        self.value = value
        self.G = 1/value
        self.Gs = None
        
        # create stamp and relative informations
        self.stamp_G = array([[1, -1], [-1, 1]])
        self.contribute = ["G"]
        self.vsource = 0

    def init_component(self, frequency):
        s = laplace(frequency)
        self.Gs = self.G / s

    def update_stamp(self, node_list, M):
        """
        Update the component's stamp in a matrix of the system's size.

        Parameters
        ----------
        node_list : list
            list of nodes.

        Returns
        -------
        None.

        """
        np = self.np
        nm = self.nm
        maxNode = npmax(node_list)
        self.stamp_G = zeros([maxNode+M, maxNode+M], dtype=complex)
        
        if np != 0:
            self.stamp_G[np - 1, np - 1] = 1
        if nm != 0:
            self.stamp_G[nm - 1, nm - 1] = 1
        if np != 0 and nm != 0:
            self.stamp_G[np - 1, nm - 1] = -1
            self.stamp_G[nm - 1, np - 1] = -1
   
    
class capacitance:
    def __init__(self, np, nm, value):
        """
        Create a capacitor component.

        Parameters
        ----------
        np : int
            positive node.
        nm : int
            negative node.
        value : float
            resustance value.

        Returns
        -------
        None.

        """
        
        self.np = np
        self.nm = nm
        self.value = value
        self.G = 1/value
        self.Gs = None
        
        # create stamp and relative informations
        self.stamp_G = array([[1, -1], [-1, 1]])
        self.contribute = ["G"]
        self.vsource = 0

    def init_component(self, frequency):
        s = laplace(frequency)
        self.Gs = self.G * s

    def update_stamp(self, node_list, M):
        """
        Update the component's stamp in a matrix of the system's size.

        Parameters
        ----------
        node_list : list
            list of nodes.

        Returns
        -------
        None.

        """
        np = self.np
        nm = self.nm
        maxNode = npmax(node_list)
        self.stamp_G = zeros([maxNode+M, maxNode+M], dtype=complex)
        
        if np != 0:
            self.stamp_G[np - 1, np - 1] = 1
        if nm != 0:
            self.stamp_G[nm - 1, nm - 1] = 1
        if np != 0 and nm != 0:
            self.stamp_G[np - 1, nm - 1] = -1
            self.stamp_G[nm - 1, np - 1] = -1