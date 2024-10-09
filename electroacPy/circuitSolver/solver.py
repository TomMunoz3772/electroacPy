#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 11:26:36 2024

@author: tom
"""
from numpy import unique, zeros, squeeze, linalg, concatenate, int16
from numpy import max as nmax
from copy import copy
from tqdm import tqdm
    
class circuit:
    def __init__(self, frequency):
        """
        Create a circuit object.

        Parameters
        ----------
        frequency : numpy array
            Range of simulation.

        Returns
        -------
        None.

        """
        self.frequency = frequency
        self.Nfft = len(frequency)
        self.G_m = None
        self.I_s = None
        self.X = None
        
        # component count
        self.node_count = 1
        self.node_list = None
        self.N = 0  # number of nodes
        self.M = 0  # number of independant voltage source
        self.Nc = 0 # total number of components
        self.components = []
        
        # id component
        self.node_id = {}
        self.block = []
        self.source_id = {}
        
    def addComponent(self, *args):
        """
        Add one or multiple component object into the circuit
    
        Parameters
        ----------
        component : component object
    
        Returns
        -------
        None.
    
        """
        
        for comp in args:
            comp.init_component(self.frequency)
            self.components.append(comp)
            self.checkNodes(comp)
                
    def checkNodes(self, component):
        """
        Looks for already existing node in circuit (stored in node_id). If none
        is found, add a new key to the node_id dictionary.

        Parameters
        ----------
        component : component object
            Electric, mechanic or acoustic components.

        Returns
        -------
        None.

        """
        
        # Check 'np' key and assign node if it doesn't exist
        if component.np != '0' and component.np not in self.node_id:
            self.node_id[component.np] = self.node_count
            self.node_count += 1
    
        # Check 'nm' key and assign node if it doesn't exist
        if component.nm != '0' and component.nm not in self.node_id:
            self.node_id[component.nm] = self.node_count
            self.node_count += 1
    
        # Check 'np1' key if it exists in the component and assign node if it doesn't exist
        if hasattr(component, 'np1') and component.np1 != '0' and component.np1 not in self.node_id:
            self.node_id[component.np1] = self.node_count
            self.node_count += 1
    
        # Check 'nm1' key if it exists in the component and assign node if it doesn't exist
        if hasattr(component, 'nm1') and component.nm1 != '0' and component.nm1 not in self.node_id:
            self.node_id[component.nm1] = self.node_count
            self.node_count += 1
                
                
    def addBlock(self, *args):
        """
        Add one or multiple block object into the circuit
    
        Parameters
        ----------
        block : block object
    
        Returns
        -------
        None.
    
        """
        for block in args:
            for comp in block.network:
                block.network[comp].init_component(self.frequency)
                self.components.append(block.network[comp])
                self.checkNodes(block.network[comp])

                
    def countNodesAndSources(self):                   
        """
        Count and identify nodes and sources present in the circuits.

        Returns
        -------
        None.

        """
        self.N = len(self.node_id)        
        self.Nc = len(self.components)
        
        for comp in self.components:
            if comp.vsource != 0:
                self.M += comp.vsource
        
        nbsource = 0
        for comp in self.components:
            if "Is" not in comp.contribute:
                comp.update_stamp(self.node_id, self.M)
            elif "Is" in comp.contribute:
                comp.update_stamp(self.node_id, self.M, nbsource)
                nbsource += comp.vsource  # add number of sources
                # print("nbsource", nbsource)
                # print("comp.np", comp.np)
                if comp.vsource == 2:
                    self.source_id[comp.np] = nbsource-1
                    self.source_id[comp.np1] = nbsource
                else:
                    self.source_id[comp.np] = nbsource
    
    def build_G(self, nf):
        G = zeros([self.N+self.M, self.M+self.N], dtype=complex)
        for comp in self.components:
            if "G" in comp.contribute:
                G += comp.stamp_G * comp.Gs[nf]
                if "Ys" in comp.__dict__:
                    G += comp.stamp_Y * comp.Ys[nf]             
        return G
    
    def build_Is(self, nf):
        I = zeros([self.N + self.M, 1], dtype=complex)
        for comp in self.components:
            if "Is" in comp.contribute:
                I += comp.stamp_I * comp.Gs[nf]
        return I
    
    def run(self):
        self.countNodesAndSources()
        self.X = zeros([self.N+self.M, self.Nfft], dtype=complex)
        
        print("Solving network...")
        for nf in tqdm(range(self.Nfft)):
            self.G_m = self.build_G(nf)
            self.I_s = self.build_Is(nf)
            self.X[:, nf] = squeeze(linalg.inv(self.G_m) @ self.I_s)
    
    def getPotential(self, node_id):
        """
        Return the potential at given node.

        Parameters
        ----------
        node_id : int or str
            Node at which the potential is probed.

        Returns
        -------
        out : numpy array
            Voltage at node "node_id".

        """

        if isinstance(node_id, int):
            node_id = str(node_id)
        node_number = self.node_id[node_id]
        
        if self.X is not None:
            if node_number <= self.N:
                out = self.X[node_number-1]
            else:
                raise Exception("Potential is given for N nodes: node_number <= N")
        return out
    
    def getFlow(self, source_id):
        """
        Return the flow at given source.

        Parameters
        ----------
        source_id : int or str
            Node of source, as defined by the user when creating the component.

        Returns
        -------
        out : numpy array
            Current at source "source_id".

        """

        if isinstance(source_id, int):
            source_id = str(source_id)
        source_number = self.source_id[source_id]
        
        if self.X is not None:
            if source_number <= self.M:
                out = self.X[self.N+source_number-1]
            else:
                raise Exception("Flow is given for M sources: source_number <= M")
    
        return out
     
    