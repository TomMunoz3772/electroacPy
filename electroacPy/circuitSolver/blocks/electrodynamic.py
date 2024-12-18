"""
Collection of electro-dynamic circuit blocks 

"""
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
        from electroacPy.circuitSolver.components.electric import resistance, inductance, capacitance
        from electroacPy.circuitSolver.components.coupler import CCVS 
        
        
        self.Re = Re
        self.Le = Le 
        self.Bl = Bl
        self.Mms = Mms
        self.Cms = Cms
        self.Rms = Rms
        self.Sd = Sd
        
        np = str(A)   # input electric
        nm = str(B)   # output electric
        np1 = str(C)  # input acoustic
        nm1 = str(D)  # output acoustic
        rnd_id = randomblock_id(3)
        
        if v_probe is not None:
            v_p = v_probe
        elif v_probe is None:
            v_p = np + "_m_8_" + rnd_id

        self.network = {"Le": inductance(np, np+"_e1_"+rnd_id, self.Le),
                        "Re": resistance(np+"_e1_"+rnd_id, np+"_e2_"+rnd_id, self.Re),
                        "Bl1": CCVS(np+"_e2_"+rnd_id, np+"_e3_"+rnd_id,
                                    np+"_m4_"+rnd_id, np+"_m5_"+rnd_id, self.Bl),
                        "Bl2": CCVS(np+"_m5_"+rnd_id, 0, 
                                    np+"_e3_"+rnd_id, nm, -self.Bl),
                        "Mms": inductance(np+"_m4_"+rnd_id, 
                                          np+"_m6_"+rnd_id, self.Mms),
                        "Cms": capacitance(np+"_m6_"+rnd_id, 
                                           np+"_m7_"+rnd_id, self.Cms),
                        "Rms": resistance(np+"_m7_"+rnd_id, 
                                          v_p, self.Rms),
                        "Sd1": CCVS(v_p, np+"_m9_"+rnd_id, 
                                    np+"_s10_"+rnd_id, np+"_s11_"+rnd_id, self.Sd),
                        "RSn": resistance(np+"_s10_"+rnd_id, np+"_s12_"+rnd_id, 1e-9),
                        "Sd2": CCVS(np+"_s12_"+rnd_id, np+"_s13_"+rnd_id, 
                                    np1, np1+"_a14_"+rnd_id, 1),
                        "Sd3": CCVS(np1+"_a14_"+rnd_id, nm1, 
                                    np+"_s13_"+rnd_id, 0, -1),
                        "Sd4": CCVS(np+"_s11_"+rnd_id, 0, 
                                    np+"_m9_"+rnd_id, 0, -self.Sd)
                        }     
        
        
class EADImport:
    def __init__(self, A, B, C, D, lpmFile, v_probe=None, rho=1.22, c=343):
        """
        Creates an electro-acoustic driver based on its input/output nodes and
        T/S parameters --- from LPM data-file.

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
        from electroacPy.circuitSolver.components.electric import resistance, inductance, capacitance
        from electroacPy.circuitSolver.components.coupler import CCVS 
        from electroacPy.speakerSim.electroAcousticDriver import loadLPM
        from numpy import arange
        
        drv = loadLPM(lpmFile, arange(10, 100, 1), c=c, rho=rho)
        
        self.Re = drv.Re
        self.Le = drv.Le 
        self.Bl = drv.Bl
        self.Mms = drv.Mms
        self.Cms = drv.Cms
        self.Rms = drv.Rms
        self.Sd = drv.Sd
        self.Vas = drv.Vas
        self.Fs = drv.Fs
        
        np = str(A)   # input electric
        nm = str(B)   # output electric
        np1 = str(C)  # input acoustic
        nm1 = str(D)  # output acoustic
        rnd_id = randomblock_id(3)
        
        if v_probe is not None:
            v_p = v_probe
        elif v_probe is None:
            v_p = np + "_m_8_" + rnd_id

        self.network = {"Le": inductance(np, np+"_e1_"+rnd_id, self.Le),
                        "Re": resistance(np+"_e1_"+rnd_id, np+"_e2_"+rnd_id, self.Re),
                        "Bl1": CCVS(np+"_e2_"+rnd_id, np+"_e3_"+rnd_id,
                                    np+"_m4_"+rnd_id, np+"_m5_"+rnd_id, self.Bl),
                        "Bl2": CCVS(np+"_m5_"+rnd_id, 0, 
                                    np+"_e3_"+rnd_id, nm, -self.Bl),
                        "Mms": inductance(np+"_m4_"+rnd_id, 
                                          np+"_m6_"+rnd_id, self.Mms),
                        "Cms": capacitance(np+"_m6_"+rnd_id, 
                                           np+"_m7_"+rnd_id, self.Cms),
                        "Rms": resistance(np+"_m7_"+rnd_id, 
                                          v_p, self.Rms),
                        "Sd1": CCVS(v_p, np+"_m9_"+rnd_id, 
                                    np+"_s10_"+rnd_id, np+"_s11_"+rnd_id, self.Sd),
                        "RSn": resistance(np+"_s10_"+rnd_id, np+"_s12_"+rnd_id, 1e-9),
                        "Sd2": CCVS(np+"_s12_"+rnd_id, np+"_s13_"+rnd_id, 
                                    np1, np1+"_a14_"+rnd_id, 1),
                        "Sd3": CCVS(np1+"_a14_"+rnd_id, nm1, 
                                    np+"_s13_"+rnd_id, 0, -1),
                        "Sd4": CCVS(np+"_s11_"+rnd_id, 0, 
                                    np+"_m9_"+rnd_id, 0, -self.Sd)
                        }     
        
        
        
        
        