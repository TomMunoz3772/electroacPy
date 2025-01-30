#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 11:32:22 2025

@author: tom
"""

import bempp.api
from bempp.api.operators.boundary import helmholtz, sparse
from bempp.api.operators.potential import helmholtz as helmholtz_potential
from bempp.api.assembly.discrete_boundary_operator import DiagonalOperator
from scipy.sparse.linalg import gmres as scipy_gmres
from bempp.api.linalg import gmres
import numpy as np
from tqdm import tqdm
import warnings
from pyopencl import CompilerWarning
import electroacPy.general as gtb

warnings.filterwarnings("ignore", message="splu requires CSC matrix format")
warnings.filterwarnings("ignore", message="splu converted its input to CSC format")
warnings.filterwarnings("ignore", category=CompilerWarning)

class pointSource:
    def __init__(self, sourcePosition, radiatingElement, velocity, frequency, 
                 c=343, rho=1.22, **kwargs):
        self.xSource = sourcePosition
        self.radiatingElement = radiatingElement
        self.velocity = velocity


