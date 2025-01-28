#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 10:40:49 2024

@author: tom
"""

from numpy import array, zeros, ones, pi, sqrt, tan, sin
from electroacPy.general.freqop import laplace


class pressureSource:
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
        
        self.np = str(np)
        self.nm = str(nm)
        self.value = value
        self.G = 1
        self.Gs = None
        
        # create stamp and relative informations
        self.stamp_G = array([[0, 0, 1], [0, 0, -1], [1, -1, 0]])
        self.stamp_I = array([[0], [0], [value]])
        self.contribute = ["G", "Is"]  
        self.vsource = 1
    
        
    def init_component(self, frequency):
        self.Gs = 1 * ones(len(frequency)) 
        
    
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
        
        maxNode = max(node_id.values())
        
        # G stamp
        self.stamp_G = zeros([maxNode+M, maxNode+M], dtype=complex)
        
        if np != 0:
            self.stamp_G[maxNode+nbsource, np-1] = 1  # sub mat B
            self.stamp_G[np-1, maxNode+nbsource] = 1  # sub mat C
        if nm != 0:
            self.stamp_G[maxNode+nbsource, nm-1] = -1 # sub mat B
            self.stamp_G[nm-1, maxNode+nbsource] = -1 # sub mat C     
        
        # I stamp
        self.stamp_I = zeros([maxNode+M, 1], dtype=complex) 
        self.stamp_I[maxNode+nbsource] = self.value #1


#%% loudspeaker relared components
class radiator:
    def __init__(self, np, nm, Sd, rho=1.22, c=343):
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
        
        self.np = str(np)
        self.nm = str(nm)
        self.Sd = Sd
        self.G = 1/Sd
        self.Gs = None
        
        self.rho = rho
        self.c = c
        self.r = sqrt(Sd/pi)
        
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
        # self.Gs = self.G * ones(len(frequency)) 
        from scipy.special import j1, struve
        
        om = 2 * pi * frequency
        k = om / self.c
        Zc = self.rho * self.c

        self.Za = (Zc*(1 - j1(2 * k * self.r) / k / self.r) +
                    1j * (Zc*struve(1, 2 * k * self.r) / k / self.r))
        self.Gs = 1/self.Za # conductance
        
        
    def update_stamp(self, node_id, M):
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
        
        maxNode = max(node_id.values())
        
        # update stamp
        self.stamp_G = zeros([maxNode+M, maxNode+M], dtype=complex)
        
        if np != 0:
            self.stamp_G[np - 1, np - 1] = 1
        if nm != 0:
            self.stamp_G[nm - 1, nm - 1] = 1
        if np != 0 and nm != 0:
            self.stamp_G[np - 1, nm - 1] = -1
            self.stamp_G[nm - 1, np - 1] = -1
    
class cavity:
    def __init__(self, np, nm, Vb, eta=1e-5, rho=1.22, c=343, losses=True):
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
        
        self.np = str(np)
        self.nm = str(nm)
        self.Vb = Vb
        self.eta = eta
        self.rho = rho
        self.c = c
        self.losses = losses
        self.G = 1
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
        s = laplace(frequency)
        # Lb = (self.Vb)**(1/3) / 3
        # Sb = self.Vb / Lb
        
        Cb = self.Vb / self.rho / self.c**2
        Rb = self.rho * self.c / self.eta / self.Vb
        if self.losses is True:
            self.Gs = s * Cb + 1/Rb # conductance
        elif self.losses is False:
            self.Gs = s * Cb
        else:
            raise Exception("Losses must be boolean.")
        
        
    def update_stamp(self, node_id, M):
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
        
        maxNode = max(node_id.values())
        
        # update stamp
        self.stamp_G = zeros([maxNode+M, maxNode+M], dtype=complex)
        
        if np != 0:
            self.stamp_G[np - 1, np - 1] = 1
        if nm != 0:
            self.stamp_G[nm - 1, nm - 1] = 1
        if np != 0 and nm != 0:
            self.stamp_G[np - 1, nm - 1] = -1
            self.stamp_G[nm - 1, np - 1] = -1

class port:
    def __init__(self, np, nm, Lp, rp, flange="single", 
                mu=1.86e-5, rho=1.22, c=343):
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
        
        self.np = str(np)
        self.nm = str(nm)
        self.Lp = Lp
        self.rp = rp
        self.Sp = pi * rp**2
        self.rho = rho
        self.c = c
        self.mu = mu
        self.G = 1
        self.Gs = None
        
        # prepare flange
        if flange == "single":
            self.flangeCoeff = 0.84
        elif flange == "both":
            self.flangeCoeff = 0.96
        elif flange == "none":
            self.flangeCoeff = 0.76
        else:
            self.flangeCoeff = 0.84
        
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
        om = 2 * pi * frequency
        s = laplace(frequency)
        
        Lt = self.Lp + self.flangeCoeff * sqrt(self.Sp)
        
        Mp = Lt*self.rho / (self.Sp)
        Rp = (2*om*self.rho*self.mu)/self.Sp * (Lt/self.rp + 1*0.7)
        Zp = s*Mp + Rp
        self.Gs = 1 / Zp # conductance
        
        
    def update_stamp(self, node_id, M):
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
        
        maxNode = max(node_id.values())
        
        # update stamp
        self.stamp_G = zeros([maxNode+M, maxNode+M], dtype=complex)
        
        if np != 0:
            self.stamp_G[np - 1, np - 1] = 1
        if nm != 0:
            self.stamp_G[nm - 1, nm - 1] = 1
        if np != 0 and nm != 0:
            self.stamp_G[np - 1, nm - 1] = -1
            self.stamp_G[nm - 1, np - 1] = -1    

class membrane:
    def __init__(self, np, nm, Cm, Mm, Rm, Sd, rho=1.22, c=343):
        self.np = str(np)
        self.nm = str(nm)
        
        self.Cm = Cm
        self.Mm = Mm
        self.Rm = Rm
        self.Sd = Sd
        self.rho = rho
        self.c = c
        self.G = 1
        self.Gs = None
        
        self.Ma = Mm / Sd**2
        self.Ca = Cm * Sd**2
        self.Ra = Rm / Sd**2
        
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
       s = laplace(frequency)
       
       Za = s*self.Ma + 1/self.Ca/s + self.Ra
       self.Gs = 1 / Za # conductance
       
       
    def update_stamp(self, node_id, M):
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
        
        maxNode = max(node_id.values())
        
        # update stamp
        self.stamp_G = zeros([maxNode+M, maxNode+M], dtype=complex)
        
        if np != 0:
            self.stamp_G[np - 1, np - 1] = 1
        if nm != 0:
            self.stamp_G[nm - 1, nm - 1] = 1
        if np != 0 and nm != 0:
            self.stamp_G[np - 1, nm - 1] = -1
            self.stamp_G[nm - 1, np - 1] = -1   
            
class closed_line:
    def __init__(self, np, nm, Lp, Sp, eta=1e-3, rho=1.22, c=343):
        """
        Add a closed cavity (nm linked to ground)

        Parameters
        ----------
        np : TYPE
            DESCRIPTION.
        nm : TYPE
            DESCRIPTION.
        Lp : TYPE
            DESCRIPTION.
        Sp : TYPE
            DESCRIPTION.
        rho : TYPE, optional
            DESCRIPTION. The default is 1.22.
        c : TYPE, optional
            DESCRIPTION. The default is 343.

        Returns
        -------
        None.

        """
        
        self.np = str(np)
        self.nm = str(nm)
        
        self.Lp = Lp
        self.Sp = Sp
        self.rp = sqrt(Sp / pi) # equivalent radius (not necessarily circular duct)
        self.Pp = 2*pi*self.rp  # equivalent perimeter (to estimate thermo-viscous losses)
        self.rho = rho
        self.c = c
        self.G = 1
        self.Gs = None
        self.eta = eta
        self.stamp_G = None
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
       om = 2 * pi * frequency
       k = om / self.c
       Zc = self.rho * self.c / self.Sp
       kl = k-1j*self.eta # kloss(self.Sp, self.Pp, k, self.rho, self.c)
       Za = Zc / 1j / tan(kl*self.Lp)
       self.Gs = 1 / Za
       
       
       
    def update_stamp(self, node_id, M):
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
        
        maxNode = max(node_id.values())
        
        # update stamp
        self.stamp_G = zeros([maxNode+M, maxNode+M], dtype=complex)
        
        if np != 0:
            self.stamp_G[np - 1, np - 1] = 1
        if nm != 0:
            self.stamp_G[nm - 1, nm - 1] = 1
        if np != 0 and nm != 0:
            self.stamp_G[np - 1, nm - 1] = -1
            self.stamp_G[nm - 1, np - 1] = -1


class open_line_T:
    def __init__(self, np0, np1, np2, nm, Lp, Sp, eta=1e-3, rho=1.22, c=343):
        self.np = str(np0) # careful with numbering and references
        self.nm = str(np1)
        self.np1 = str(np2)
        self.nm1 = str(nm)  # ground
        self.Lp = Lp
        self.Sp = Sp  # surface
        self.rp = sqrt(Sp / pi)  # equivalent radius (not necessarily circular duct)
        self.Pp = 2 * pi * self.rp  # equivalent perimeter (to estimate thermo-viscous losses)
        self.rho = rho
        self.c = c   
        self.eta = eta
        
        self.G = 1
        self.Gs = None
        self.Ys = None
        self.stamp_G = None
        self.stamp_Y = None
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
       
       LL = self.Lp #+ 0.73*self.rp*2
       om = 2 * pi * frequency
       k = om / self.c
       Zc = self.rho * self.c / self.Sp
       kl = k - 1j*self.eta #kloss(self.Sp, self.Pp, k, self.rho, self.c) # k - 1j*self.eta #
       Zt =  1j * Zc * tan(kl*LL/2)
       Yt = 1j * sin(kl*LL) / Zc
       self.Gs = 1 / Zt
       self.Ys = Yt

    def update_stamp(self, node_id, M):
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
            np0 = node_id[self.np]
        else:
            np0 = 0
            
        if self.nm != '0':
            np1 = node_id[self.nm]
        else:
            np1 = 0
        
        if self.np1 != '0':
            np2 = node_id[self.np1]
        else:
            np2 = 0
            
        if self.nm1 != '0':
            nm = node_id[self.nm1]
        else:
            nm = 0
        
        maxNode = max(node_id.values())                

        self.stamp_G = zeros([maxNode+M, maxNode+M], dtype=complex)
        self.stamp_Y = zeros([maxNode+M, maxNode+M], dtype=complex)
        
        if np0 != 0:
            self.stamp_G[np0 - 1, np0 - 1] = 1
        if np1 != 0:
            self.stamp_G[np1 - 1, np1 - 1] = 2
            self.stamp_Y[np1 - 1, np1 - 1] = 1
        if np2 != 0:
            self.stamp_G[np2 - 1, np2 - 1] = 1
        if nm != 0:
            self.stamp_Y[nm - 1, nm - 1] = 1
            
        if np0 != 0 and np1 != 0:
            self.stamp_G[np0 - 1, np1 - 1] = -1
            self.stamp_G[np1 - 1, np0 - 1] = -1 
        
        if np1 != 0 and np2 != 0:
            self.stamp_G[np1 - 1, np2 - 1] = -1
            self.stamp_G[np2 - 1, np1 - 1] = -1 
            
        if np1 != 0 and nm !=0:
            self.stamp_Y[np1 - 1, nm - 1] = -1
            self.stamp_Y[nm - 1, np1 - 1] = -1 


def kloss(section, perimeter, k, rho=1.22, c=343):
    """
    Losses as expressed by Jean-Pierre Dalmont in his university courses

    Input parameters
    -----
    section: duct cross-sectional area
    perimeter: perimeter of the cross-section
    k: input wavenumber

    ------
    return: k_l --> wavenumber with losses
    """
    f = k * c / 2 / pi
    om = 2 * pi * f
    alpha = sqrt(f) * (0.95e-5 + 2.03e-5) * perimeter / 2 / section  # +1 was added
    k_l = om / c * (1 + (1 - 1j) * alpha)
    return k_l


