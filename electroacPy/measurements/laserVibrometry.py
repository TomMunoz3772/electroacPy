import numpy as np
from electroacPy.general import freqop, geometry, io
import matplotlib.pyplot as plt


class laserVibrometry_UFF:
    def __init__(self, file_path, spatial_rotation, freq_array, 
                 useAverage=False, inputVoltage=1):
        
        """
        Import UFF file to load acceleration data.
        """
        
        self.file_path = file_path
        self.spatial_rotation = spatial_rotation
        self.identifier = "PLV"

        # load and decimate data to match frequency axis
        freq_meas, acc_meas, v_meas, point_cloud = io.loadUFF(file_path)
        freq_dec, index_dec = freqop.decimate_frequency_axis(freq_meas, freq_array)
        acc_dec = acc_meas[:, index_dec] * inputVoltage
        v_dec = v_meas[:, index_dec] * inputVoltage

        self.acc_dec = acc_dec
        self.freq_dec = freq_dec
        
        # useful data
        if np.all(v_dec == 0) is True or np.array(v_dec).size == 0:
            self.v_point = acc_dec / freqop.laplace(freq_dec)       # velocity at each points (nPoints, Nfft)
        else:
            self.v_point = acc_dec / freqop.laplace(freq_dec) # v_dec
        self.v           = np.ones(len(freq_array), dtype=complex)  # decoy for post processing ?
        self.point_cloud = geometry.rotatePointCloud(geometry.recenterZero(point_cloud), spatial_rotation)

        self.v_mean = np.zeros([len(freq_dec)], dtype=complex)
        for f in range(len(freq_dec)):
            self.v_mean[f] = np.mean(self.v_point[:, f])

        if useAverage is True:
            for f in range(len(freq_dec)):
                self.v_point[:, f] = self.v_mean[f]

        # identifiers
        self.ref2bem   = None
        self.poly_data = True  # is class from polytech?

    def plot_pointCloud(self):
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        ax.scatter(self.point_cloud[1:, 0]*1e2, self.point_cloud[1:, 1]*1e2, self.point_cloud[1:, 2]*1e2)
        ax.scatter(self.point_cloud[0, 0]*1e2, self.point_cloud[0, 1]*1e2, self.point_cloud[0, 2]*1e2, 'r')
        ax.axis('equal')
        ax.set(xlabel='x [cm]', ylabel='y [cm]', zlabel='z [cm]', title='Vibrometry points')
        return plt.show()

    def plot_velocity(self, point_index):
        fig, ax = plt.subplots(2, 1)
        ax[0].semilogx(self.freq_dec, np.abs(self.v_point[point_index, :]))
        ax[1].semilogx(self.freq_dec, np.angle(self.v_point[point_index, :]))
        for i in range(2):
            ax[i].set(xlabel='Frequency [Hz]', ylabel='Velocity [m/s]')
            ax[i].grid(which='both')
        plt.tight_layout()
        return plt.show()

    def plotXVA(self):
        """
        Plot the displacement, velocity, and acceleration frequency responses.

        Returns
        -------
        None

        """
        s = freqop.laplace(self.freq_dec)
        Hx = self.v_mean / s
        Hv = self.v_mean
        Ha = self.v_mean * s

        fig, ax = plt.subplots(3, 1)
        ax[0].semilogx(self.freq_dec, np.abs(Hx * 1e3), label='Displacement')
        ax[1].semilogx(self.freq_dec, np.abs(Hv), label='Velocity')
        ax[2].semilogx(self.freq_dec, np.abs(Ha), label='Acceleration')
        ax[2].set(xlabel="Frequency [Hz]")
        ax[0].set(ylabel="mm")
        ax[1].set(ylabel="m/s")
        ax[2].set(ylabel="m/s^2")
        for i in range(3):
            ax[i].grid(which='both')
            ax[i].legend(loc='best')
        plt.tight_layout()
        return plt.show()



class laserVibrometry:
    def __init__(self, Hv, X, inputVoltage=1):
        """
        
        Import user-defined vib data in numpy arrays.        
        
        Parameters
        ----------
        Hv : numpy array
            Velocity data (nPoints, NFFT). Nfft should corresponds to the number of points
            in your simulation.
        X : numpy array
            Position data (nPoints, 3) in cartesian coordinates.
    
        Returns
        -------
        None.

        """
        self.identifier = "PLV"
        self.v_point = Hv * inputVoltage 
        self.point_cloud = X
        
        self.v = np.ones(Hv.shape[1], dtype=complex) 
        self.v_mean = np.mean(self.v_point, 1)

        # identifiers
        self.ref2bem   = None
        self.poly_data = True  # is class from polytech?

    def plot_pointCloud(self):
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        ax.scatter(self.point_cloud[1:, 0]*1e2, self.point_cloud[1:, 1]*1e2, self.point_cloud[1:, 2]*1e2)
        ax.scatter(self.point_cloud[0, 0]*1e2, self.point_cloud[0, 1]*1e2, self.point_cloud[0, 2]*1e2, 'r')
        ax.axis('equal')
        ax.set(xlabel='x [cm]', ylabel='y [cm]', zlabel='z [cm]', title='Vibrometry points')
        return plt.show()

    def plot_velocity(self, point_index):
        fig, ax = plt.subplots(2, 1)
        ax[0].semilogx(self.freq_dec, np.abs(self.v_point[point_index, :]))
        ax[1].semilogx(self.freq_dec, np.angle(self.v_point[point_index, :]))
        for i in range(2):
            ax[i].set(xlabel='Frequency [Hz]', ylabel='Velocity [m/s]')
            ax[i].grid(which='both')
        plt.tight_layout()
        return plt.show()