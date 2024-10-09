"""
Collection of acoustic circuit blocks 

"""

from electroacPy.circuitSolver.components.electric import resistance, capacitance, inductance
from electroacPy.circuitSolver.components.acoustic import radiator
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
    def __init__(self, A, Vb, Lp, Sp, fs, Vas, Cas,
                 Ql=10, Qa=100, k=0.7, rho=1.22, c=343, p_probe=None):
        """
        Create a ported enclosure. Uses loudspeaker compliance to determine 
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
        k: float,
            length correction. (one flanged termination) 0.6 < k < 0.9 (both termination are flanged)

        Returns
        -------
        None.kwargs
        """
        
        np = str(A)
        nm = 0
        rnd_id = randomblock_id(3)

        # parameter computation
        Cab = Vb / rho / c**2
        Lt  = Lp + k * sqrt(Sp/pi)
        Mp  = rho*Lt/Sp
        
        fb = 1 / (2*pi*sqrt(Mp*Cab))
        wb = 2*pi*fb
        
        Rab = 1 / wb / Qa / Cab
        Ral = Ql / wb / Cab
        Rp  = 1 / wb / 100 / Cab
        
        if p_probe is not None:
            np_rad = p_probe
        else:
            np_rad = np+"_3_"+rnd_id
        
        self.network = {"Rab": resistance(np, np+"_1_"+rnd_id, Rab),
                        "Cab": capacitance(np+"_1_"+rnd_id, nm, Cab),
                        "Ral": resistance(np, nm, Ral),
                        "Mp" : inductance(np, np+"_2_"+rnd_id, Mp),
                        "Rp" : resistance(np+"_2_"+rnd_id, np_rad, Rp),
                        "Rad": radiator(np_rad, nm, Sp, rho=rho, c=c),
                        }