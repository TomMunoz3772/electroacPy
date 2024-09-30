"""
Collection of electro-dynamic circuit blocks 

"""

from electroacPy.circuitSolver.components.electric import resistance, inductance, capacitance
from electroacPy.circuitSolver.components.coupler import CCVS 
import random, string

def randomblock_id(length):
   letters = string.ascii_lowercase
   return ''.join(random.choice(letters) for i in range(length))

class EAD:
    def __init__(self, A, B, C, D, 
                 Le, Re, Cms, Mms, Rms, Bl, Sd, v_probe=None):
        """
        Creates an electro-acoustic driver based on its input/output nodes and
        T/S parameters.

        Parameters
        ----------
        A : int or str
            positive electrical node.
        B : int or str
            negative electrical node..
        C : int or str
            positive acoustical node..
        D : int or str
            negative acoustical node..
        Le : float
            coil inductance.
        Re : float
            coil's electrical resistance.
        Cms : float
            suspensions compliance.
        Mms : float
            moving mass.
        Rms : float
            mechanical losses.
        Bl : float
            electrical to mechanical force factor.
        Sd : float
            effective radiating surface.

        Returns
        -------
        None.

        """
        
        self.Re = Re
        self.Le = Le 
        self.Bl = Bl
        self.Mms = Mms
        self.Cms = Cms
        self.Rms = Rms
        self.Sd = Sd
        
        np = str(A)   # input electric
        nm = str(B) # input acoustic
        np1 = str(C)   # output electric
        nm1 = str(D) # output acoustic
        rnd_id = randomblock_id(3)
        
        if v_probe is not None:
            v_p = v_probe
        elif v_probe is None:
            v_p = "m_4" + rnd_id
        
        
        # self.network = {"Le": inductance(np, np+"_e1_"+rnd_id, self.Le),
        #                 "Re": resistance(np+"_e1_"+rnd_id, np+"_e2_"+rnd_id, self.Re),
        #                 "Bl1": CCVS(np+"_e2_"+rnd_id, np+"_e3_"+rnd_id,
        #                             np+"_m4_"+rnd_id, np+"_m5_"+rnd_id, self.Bl),
        #                 "Bl2": CCVS(np+"_m5_"+rnd_id, 0, 
        #                             nm, np+"_e3_"+rnd_id, self.Bl),
        #                 "Mms": inductance(np+"_m4_"+rnd_id, 
        #                                   np+"_m6_"+rnd_id, self.Mms),
        #                 "Cms": capacitance(np+"_m6_"+rnd_id, 
        #                                    np+"_m7_"+rnd_id, self.Cms),
        #                 "Rms": resistance(np+"_m7_"+rnd_id, 
        #                                   np+"_m8_"+rnd_id, self.Rms),
        #                 "Sd1": CCVS(np+"_m8_"+rnd_id, np+"_m9_"+rnd_id, 
        #                             np+"_s10_"+rnd_id, np+"_s11_"+rnd_id, self.Sd),
        #                 "Sd2": CCVS(np+"_s10_"+rnd_id, np+"_s12_"+rnd_id, 
        #                             np1, np1+"_a14_"+rnd_id, 1),
        #                 "Sd3": CCVS(np1+"_a14_"+rnd_id, nm1, 0, 
        #                             np+"_s12_"+rnd_id, 1),
        #                 "Sd4": CCVS(np+"_s11_"+rnd_id, 0, 0, 
        #                             np+"_m9_"+rnd_id, self.Sd)
        #                 }
        
        self.network = {"Le": inductance(np, np+"_e1_"+rnd_id, self.Le),
                        "Re": resistance(np+"_e1_"+rnd_id, np+"_e2_"+rnd_id, self.Re),
                        "Bl1": CCVS(np+"_e2_"+rnd_id, np+"_e3_"+rnd_id,
                                    np+v_p, np+"_m5_"+rnd_id, self.Bl),
                        "Bl2": CCVS(np+"_m5_"+rnd_id, 0, 
                                    nm, np+"_e3_"+rnd_id, self.Bl),
                        "Mms": inductance(np+v_p, 
                                          np+"_m6_"+rnd_id, self.Mms),
                        "Cms": capacitance(np+"_m6_"+rnd_id, 
                                           np+"_m7_"+rnd_id, self.Cms),
                        "Rms": resistance(np+"_m7_"+rnd_id, 
                                          np+"_m8_"+rnd_id, self.Rms),
                        "Sd1": CCVS(np+"_m8_"+rnd_id, np+"_m9_"+rnd_id, 
                                    np+"_s10_"+rnd_id, np+"_s11_"+rnd_id, self.Sd),
                        "Sd2": CCVS(np+"_s10_"+rnd_id, np+"_s12_"+rnd_id, 
                                    np1, np1+"_a14_"+rnd_id, 1),
                        "Sd3": CCVS(np1+"_a14_"+rnd_id, nm1, 0, 
                                    np+"_s12_"+rnd_id, 1),
                        "Sd4": CCVS(np+"_s11_"+rnd_id, 0, 0, 
                                    np+"_m9_"+rnd_id, self.Sd)
                        }
