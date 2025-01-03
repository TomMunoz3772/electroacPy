#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 10:40:49 2024

@author: tom
"""

from numpy import array, zeros, ones

class CCVS:
    def __init__(self, np, nm, np1, nm1, value):
        """
        Create a gyrator component

        Parameters
        ----------
        np : TYPE
            DESCRIPTION.
        nm : TYPE
            DESCRIPTION.
        np1 : TYPE
            DESCRIPTION.
        nm1 : TYPE
            DESCRIPTION.
        value : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        
        self.np = str(np)
        self.nm = str(nm)
        self.np1 = str(np1)
        self.nm1 = str(nm1)
        self.value = value
        self.G = 1
        self.Gs = None
        
        # create stamp and relative informations
        self.stamp_G = array([[0, 0, 1], 
                              [0, 0, -1], 
                              [1, -1, 0]])
        self.stamp_I = None
        self.contribute = ["G", "Is"]  
        self.vsource = 2
        
    
    def init_component(self, frequency):
        self.Gs = ones(len(frequency))
    
    def update_stamp(self, node_id, M, nbsource):
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
        if self.np != '0':
            np = node_id[self.np]
        else:
            np = 0
            
        if self.nm != '0':
            nm = node_id[self.nm]
        else:
            nm = 0
        
        if self.np1 != '0':
            np1 = node_id[self.np1]
        else:
            np1 = 0
            
        if self.nm1 != '0':
            nm1 = node_id[self.nm1]
        else:
            nm1 = 0
        
        maxNode = max(node_id.values())

        # G stamp
        self.stamp_G = zeros([maxNode+M, maxNode+M], dtype=complex)
        
        # current sensor
        if np != 0:
            self.stamp_G[np-1, maxNode+nbsource] = 1
            self.stamp_G[maxNode+nbsource, np-1] = 1
        if nm != 0:
            self.stamp_G[nm-1, maxNode+nbsource] = -1
            self.stamp_G[maxNode+nbsource, nm-1] = -1        
        
        # secondary source
        if np1 != 0:
            self.stamp_G[np1-1, maxNode+nbsource+1] = 1
            self.stamp_G[maxNode+nbsource+1, np1-1] = 1
        if nm1 != 0:
            self.stamp_G[nm1-1, maxNode+nbsource+1] = -1
            self.stamp_G[maxNode+nbsource+1, nm1-1] = -1 
        
        # Bl factor
        self.stamp_G[maxNode+nbsource+1, maxNode+nbsource] = -self.value
                
        # Y stamp
        self.stamp_I = zeros([maxNode+M, 1], dtype=complex)
        
        