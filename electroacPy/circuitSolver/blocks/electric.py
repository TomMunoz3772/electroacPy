"""
Collection of electric circuit blocks

"""

from electroacPy.circuitSolver.components.electric import inductance, capacitance
import random, string
from numpy import pi, sqrt

def randomblock_id(length):
   letters = string.ascii_lowercase
   return ''.join(random.choice(letters) for i in range(length))


class lowpass_butter:
    def __init__(self, A, B, order, fc, Re):
        """
        Create a lowpass-filter block. From Leo Beranek and Tim Mellow - 
        Acoustics Sound Fields, Transducers and Vibration.
        
        Parameters
        ----------
        np : int or str
            positive connection.
        nm : int or str
            negative connection.
        order : int
            order of the low-pass filter.
        fc : float
            cut-off frequency - up to order 6.
        Re : float
            Resistance of connected loudspeaker 
            (taken into account for component selection).

        Returns
        -------
        Values of capacitance and inductance.

        """
        self.A = A
        self.B = B
        self.order = order
        self.fc = fc
        self.Re = Re
        
        om0 = 2 * pi *self.fc
        L = Re / om0
        C = 1 / om0 / Re
        rnd_id = randomblock_id(3)
        np = str(A)
        nm = str(B)
        
        self.network = None

        if self.order == 1:
            self.network = {"L1": inductance(np, nm, L)}
            self.L = L
        
        elif self.order == 2:
            L1 = 2 * L
            C1 = 0.5 * C
            self.network = {"L1": inductance(np, nm, 2*L),
                            "C1": capacitance(nm, 0, 0.5*C)}
            self.L = L1
            self.C = C1 
        
        elif self.order == 3:
            L1 = 3/2 * L
            L2 = 1/2 * L
            C1 = 4/3 * C
            self.network = {"L1": inductance(np, np+'_lp3_'+rnd_id, 3/2*L),
                            "C1": capacitance(np+'_lp3_'+rnd_id, 0, 4/3*C),
                            "L2": inductance(np+'_lp3_'+rnd_id, nm, 1/2*L)}
            self.L = [L1, L2]
            self.C = [C1]
        
        elif self.order == 4:
            L1 = 4*sqrt(2)/3*L
            L2 = 2*sqrt(2)/3*L
            C1 = 9/4/sqrt(2)*C
            C2 = 1/2/sqrt(2)*C
            self.network = {"L1": inductance(np, np+'_lp3_'+rnd_id, L1),
                            "C1": capacitance(np+'_lp3_'+rnd_id, 0, C1),
                            "L2": inductance(np+'_lp3_'+rnd_id, nm, L2),
                            "C2": capacitance(nm, 0, C2)}
            self.L = [L1, L2]
            self.C = [C1, C2]
        
        elif self.order == 5:
            L1 = 5/(sqrt(5)+1)*L
            L2 = 2 * sqrt(5) / (sqrt(5) + 1) * L
            L3 = 1 / (sqrt(5) + 1) * L
            C1 = (2/sqrt(5) + 1) * 2 / sqrt(5) * C
            C2 = 2 / sqrt(5) * C
            self.network = {"L1": inductance(np, np+'_1lp3_'+rnd_id, L1),
                            "C1": capacitance(np+'_1lp3_'+rnd_id, 0, C1),
                            "L2": inductance(np+'_1lp3_'+rnd_id, 
                                              np+'_2lp3_'+rnd_id, L2),
                            "C2": capacitance(np+'_2lp3_'+rnd_id, 0, C2),
                            "L3": inductance(np+'_2lp3_'+rnd_id, nm, L3)}
            self.L = [L1, L2, L3]
            self.C = [C1, C2]
        
        elif self.order == 6:
            L1 = 9/5*L
            L2 = 81/55 * L
            L3 = 8/11 * L
            C1 = 50/27 * C
            C2 = 121/108 * C
            C3 = 1/4 * C
            self.network = {"L1": inductance(np, np+'_1lp3_'+rnd_id, L1),
                            "C1": capacitance(np+'_1lp3_'+rnd_id, 0, C1),
                            "L2": inductance(np+'_1lp3_'+rnd_id, 
                                              np+'_2lp3_'+rnd_id, L2),
                            "C2": capacitance(np+'_2lp3_'+rnd_id, 0, C2),
                            "L3": inductance(np+'_2lp3_'+rnd_id, nm, L3),
                            "C3": capacitance(nm, 0, C3)}
            self.L = [L1, L2, L3]
            self.C = [C1, C2, C3]
        else:
            raise Exception("Order cannot be higher than 6")

class highpass_butter:
    def __init__(self, A, B, order, fc, Re):
        """
        Create a highpass-filter block. From Beranek and Tim Mellow - 
        Acoustics Sound Fields, Transducers and Vibration.

        Parameters
        ----------
        np : int or str
            positive connection.
        nm : int or str
            negative connection.
        order : int
            order of the low-pass filter.
        fc : float
            cut-off frequency - up to order 6.
        Re : float
            Resistance of connected loudspeaker 
            (taken into account for component selection).

        Returns
        -------
        Values of capacitance and inductance.

        """
        self.A = A
        self.B = B
        self.order = order
        self.fc = fc
        self.Re = Re
        
        om0 = 2 * pi *self.fc
        L = Re / om0
        C = 1 / om0 / Re
        rnd_id = randomblock_id(3)
        np = str(A)
        nm = str(B)

        if self.order == 1:
            self.network = {"C1": capacitance(np, nm, C)}
            self.C = C
        
        elif self.order == 2:
            L1 = 2 * L
            C1 = 0.5 * C
            self.network = {"C1": capacitance(np, nm, C1),
                            "L1": inductance(nm, 0, L1)}
            self.L = L1
            self.C = C1 
        
        elif self.order == 3:
            C1 = 2/3 * C
            C2 = 2 * C
            L1 = 3/4 * L
            self.network = {"C1": capacitance(np, np+'_lp3_'+rnd_id, C1),
                            "L1": inductance(np+'_lp3_'+rnd_id, 0, L1),
                            "C2": capacitance(np+'_lp3_'+rnd_id, nm, C2)}
            self.L = [L1]
            self.C = [C1, C2]
        
        elif self.order == 4:
            C1 = 3/4/sqrt(2)*C
            C2 = 3/2/sqrt(2)*C
            L1 = 4*sqrt(2)/9*L
            L2 = 2*sqrt(2)*L
            self.network = {"C1": capacitance(np, np+'_lp3_'+rnd_id, C1),
                            "L1": inductance(np+'_lp3_'+rnd_id, 0, L1),
                            "C2": capacitance(np+'_lp3_'+rnd_id, nm, C2),
                            "L2": inductance(nm, 0, L2)}
            self.L = [L1, L2]
            self.C = [C1, C2]
        
        elif self.order == 5:
            C1 = (sqrt(5)+1)/5 * C
            C2 = (sqrt(5)+1)/2/sqrt(5) * C
            C3 = (sqrt(5)+1) * C
            L1 = 1/(2/sqrt(5)+1) * sqrt(5)/2 * L
            L2 = sqrt(5)/2 * L
            self.network = {"C1": capacitance(np, np+'_1lp3_'+rnd_id, C1),
                            "L1": inductance(np+'_1lp3_'+rnd_id, 0, L1),
                            "C2": capacitance(np+'_1lp3_'+rnd_id, 
                                              np+'_2lp3_'+rnd_id, C2),
                            "L2": inductance(np+'_2lp3_'+rnd_id, 0, L2),
                            "C3": capacitance(np+'_2lp3_'+rnd_id, nm, C3)}
            self.L = [L1, L2]
            self.C = [C1, C2, C3]
        
        elif self.order == 6:
            C1 = 9/5 * C
            C2 = 55/81 * C
            C3 = 11/8 * C
            L1 = 27/50*L
            L2 = 108/121 * L
            L3 = 4 * L
            self.network = {"C1": capacitance(np, np+'_1lp3_'+rnd_id, C1),
                            "L1": inductance(np+'_1lp3_'+rnd_id, 0, L1),
                            "C2": capacitance(np+'_1lp3_'+rnd_id, 
                                              np+'_2lp3_'+rnd_id, C2),
                            "L2": inductance(np+'_2lp3_'+rnd_id, 0, L2),
                            "C3": capacitance(np+'_2lp3_'+rnd_id, nm, C3),
                            "L3": inductance(nm, 0, L3)}
            self.L = [L1, L2, L3]
            self.C = [C1, C2, C3]
        else:
            raise Exception("Order cannot be higher than 6.")