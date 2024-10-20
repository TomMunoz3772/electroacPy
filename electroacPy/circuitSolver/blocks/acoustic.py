"""
Collection of acoustic circuit blocks 

"""

from electroacPy.circuitSolver.components.electric import resistance, capacitance, inductance
from electroacPy.circuitSolver.components.acoustic import cavity, port, open_line_T
from numpy import sqrt, pi
import random, string



def randomblock_id(length):
   letters = string.ascii_lowercase
   return ''.join(random.choice(letters) for i in range(length))

class sealedEnclosure:
    def __init__(self, A, Vb, fs, Vas, Cas,
                 Ql=10, Qa=100, rho=1.22, c=343):
        """
        Create a sealed enclosure. Uses loudspeaker compliance to determine 
        losses.

        Parameters
        ----------
        A : int or str,
            Input connection.
        Vb : float,
            enclosure volume.
        fs: float, 
            resonance frequency of drive unit.
        Vas: float, 
            equivalent volume of drive unit.
        Cas: float, 
            compliance of suspension in acoustic domain.
        Ql: float,
            Q factor related to leaks. (low leaks) 5 < Ql < 30 (high leaks)
        Qa: float,
            Q factor related to damping in the enclosure. (high-damping) 5 < Qa < +100 (low damping)
            if Q = 1/Cab -> no damping at all.
        
    

        Returns
        -------
        None.kwargs

        """
                
        np = str(A)
        nm = 0
        rnd_id = randomblock_id(3)

        # parameter computation
        fb = fs * sqrt(1 + Vas / Vb)
        wb = 2*pi*fb
        
        Cab = Vb / rho / c**2
        Rab = 1 / wb / Qa  / Cab
        Ral = Ql / wb / Cab
        
        self.network = {"Rab": resistance(np, np+rnd_id, Rab),
                        "Cab": capacitance(np+rnd_id, nm, Cab),
                        "Ral": resistance(np, nm, Ral)}
        
        
class portedEnclosure:
    def __init__(self, A, B, Vb, Lp, rp,
                 Ql=10, mu=1.86e-5, rho=1.22, c=343):
        """
        Create a ported enclosure. Uses loudspeaker compliance to determine 
        losses.

        Parameters
        ----------
        A : int or str,
            Input connection.
        Vb : float,
            enclosure volume.
        Ql: float,
            Q factor related to leaks. (low leaks) 5 < Ql < 30 (high leaks)
        Qa: float,
            Q factor related to damping in the enclosure. (high-damping) 5 < Qa < +100 (low damping)
            if Q = 1/Cab -> no damping at all.
        k: float,
            length correction. (one flanged termination) 0.6 < k < 0.9 (both termination are flanged)

        Returns
        -------
        None.kwargs
        """
        
        np = str(A)
        nm = str(B)
        rnd_id = randomblock_id(3)

        Sp = pi*rp**2
        
        # parameter computation
        Cab = Vb / rho / c**2
        Lt  = Lp + 0.73*rp*2     # Lp + k * sqrt(Sp/pi)
        Mp  = Lt * rho / Sp
        
        fb = 1 / (2*pi*sqrt(Mp*Cab))
        wb = 2*pi*fb
        
        Ral = Ql / wb / Cab
        # Ral = rho * c / 1e-5 / Vb

        print("Ral = ", Ral)
        print("Mp = ", Mp)
        print("Rp: ",  (2*10*2*pi*rho*mu)/Sp * (Lt/rp + 1*0.7))
        
        self.network = {"Cab": capacitance(np, 0, Cab),
                        "Ral": resistance(np, 0, Ral),
                        "Port": port(np, nm, Lp, rp, rho, c, mu)}
                        # "Port": open_line_T(np, np+"_1_"+rnd_id, nm, 0, 
                        #                     Lt, Sp)}
        
    