import numpy as np
from matplotlib import pyplot as plt
import generalToolbox as gtb
import electroacPy
from tqdm import tqdm

## import fenicsx modules and submodules
import ufl
from dolfinx import geometry
from dolfinx.fem import Function, FunctionSpace, assemble_scalar, form, Constant, dirichletbc
from dolfinx.fem.petsc import LinearProblem
from dolfinx.io import XDMFFile
from dolfinx.mesh import create_unit_square
from ufl import dx, grad, inner
from mpi4py import MPI
from petsc4py import PETSc
import pyvista
from dolfinx import plot
from dolfinx.io import gmshio


class fem:
    def __init__(self, mshPath,
                 radiatingSurfaces,
                 surfaceVelocity,
                 freq_array,
                 xMic=None,
                 solver_options=None,
                 c=343,
                 rho=1.22):

        self.mshPath = mshPath
        self.radiatingSurfaces = radiatingSurfaces
        self.surfaceVelocity = surfaceVelocity
        self.freq_array = freq_array
        self.xMic = np.array(xMic)
        self.rho = rho
        self.c = c
        self.solver_options = solver_options

        # if xMic == None:
        #     self.nMic = 0
        # else:
        #     self.nMic = len(xMic)
        #     if self.nMic == 1:
        #         self.xMic = np.reshape(self.xMic, [1, 3])

        ## definition of the study
        self.msh, self.cell_tags, self.facet_tags = gmshio.read_from_msh(mshPath, MPI.COMM_SELF, 0, gdim=3)
        self.V = FunctionSpace(self.msh, ('Lagrange', 1))
        u = ufl.TrialFunction(self.V)  # test function
        v = ufl.TestFunction(self.V)  # trial function
        f = Function(self.V)

        self.ds = ufl.Measure('ds', domain=self.msh, subdomain_data=self.facet_tags)  # integration over radiating surfaces

        # The flux is 0 for all rigid walls
        # g_rig = Constant(msh, PETSc.ScalarType(0))

        # Definition of the variables of the problem
        frequency = 100
        omega = 2 * np.pi * frequency
        k = omega / c
        self.vn = 10 + 10j  # normal velocity

        # equation (weak form stuff)
        self.k_sq = Constant(self.msh, PETSc.ScalarType(k ** 2))
        a = (ufl.inner(ufl.nabla_grad(u), ufl.nabla_grad(v)) * ufl.dx - self.k_sq * ufl.inner(u, v) * ufl.dx)

        # to update
        self.g_in = Constant(self.msh, 1j * omega * rho * 1)
        self.radSurf_tmp = Constant(self.msh, PETSc.ScalarType(self.radiatingSurfaces[0]))



        L = ufl.inner(self.g_in, v) * self.ds(int(self.radSurf_tmp.value.real))

        self.uh = Function(self.V)  # initialization of the solution
        if solver_options != None:
            self.problem = LinearProblem(a, L, u=self.uh, petsc_options=self.solver_options)
            # self.problem.solve()
        else:
            self.problem = LinearProblem(a, L, u=self.uh)
            # self.problem.solve()

        # radiating surfaces (m^2) and solution output (array creations)
        self.p_surf = np.zeros([len(radiatingSurfaces), len(freq_array)],
                               dtype=complex)  # pressure on radiating surfaces
        self.S = np.zeros(len(radiatingSurfaces), dtype=complex)
        # self.Zin = np.zeros([len(radiatingSurfaces), len(freq_array)], dtype=complex)
        for i in range(len(radiatingSurfaces)):
            self.S[i] = assemble_scalar(form(1 * self.ds(radiatingSurfaces[i])))

        self.solution_mat = np.zeros((len(radiatingSurfaces), len(freq_array), len(self.uh.x.array), 3))
        # self.p_mic_mat = np.zeros((len(self.radiatingSurfaces), len(self.freq_array), self.nMic), dtype=complex)
        # self.p_mic = np.zeros((len(self.freq_array), self.nMic), dtype=complex)

    def solve(self):
        print("Solving pressure")
        for nf in tqdm(range(0, len(self.freq_array))):
            freq = self.freq_array[nf]
            omega = freq * 2 * np.pi

            for nspk in range(len(self.radiatingSurfaces)): # run through all radiating surfaces
                self.g_in.value = 1j * omega * self.rho * self.surfaceVelocity[nspk, nf]
                self.radSurf_tmp.value = self.radiatingSurfaces[nspk]

                # k with losses
                k = omega / self.c
                self.k_sq.value = k ** 2
                self.problem.solve()

                # ---------------------------------------------------------------
                # Microphone pressure at specified point evaluation
                # if self.nMic != 0:
                #     for mic in range(self.nMic):
                #         points = np.zeros((3, 1))
                #         points[0][0] = self.xMic[mic, 0]
                #         points[1][0] = self.xMic[mic, 1]
                #         points[2][0] = self.xMic[mic, 2]
                #
                #         bb_tree = geometry.BoundingBoxTree(self.msh, self.msh.topology.dim)
                #         cells = []
                #         points_on_proc = []
                #         # Find cells whose bounding-box collide with the points
                #         cell_candidates = geometry.compute_collisions(bb_tree, points.T)
                #         # Choose one of the cells that contains the point
                #         colliding_cells = geometry.compute_colliding_cells(self.msh, cell_candidates, points.T)
                #
                #         for i, point in enumerate(points.T):
                #             if len(colliding_cells.links(i)) > 0:
                #                 points_on_proc.append(point)
                #                 cells.append(colliding_cells.links(i)[0])
                #         points_on_proc = np.array(points_on_proc, dtype=np.float64)
                #         u_values = self.uh.eval(points_on_proc, cells)
                #
                #         # building of sound pressure vector
                #         self.p_mic_mat[nspk, nf, mic] = u_values
                # ---------------------------------------------------------------


                # input impedance
                # for i in range(len(self.radiatingSurfaces)):
                #     self.p_surf[i, nf] = assemble_scalar(form(self.uh * self.ds(self.radiatingSurfaces[i]))) / self.S[i]
                #     Q = self.S[i] * self.vn
                #     self.Zin[i, nf] = self.p_surf[i, nf] / Q

                self.solution_mat[nspk, nf, :, 0] = self.uh.x.array.real
                self.solution_mat[nspk, nf, :, 1] = self.uh.x.array.imag

            # reshape output
            self.solution = np.zeros((len(self.freq_array), len(self.uh.x.array), 3))
            for spk in range(len(self.radiatingSurfaces)):
                self.solution += self.solution_mat[nspk, :, :, :]
                # if type(self.xMic) == np.ndarray:
                #     self.p_mic += self.p_mic_mat[nspk, :, :]
        return None

    def plotSolution(self, frequency, output='imag'):
        idxF, Fvalue = gtb.findInArray(self.freq_array, frequency)
        u_topology, u_cell_types, u_geometry = plot.create_vtk_mesh(self.V)

        u_grid = pyvista.UnstructuredGrid(u_topology, u_cell_types, u_geometry)

        if output == "imag":
            strData = "data p imag"
            plot_cmap = 'seismic'
            u_grid.point_data[strData] = self.solution[idxF, :, 1]
        elif output == 'real':
            strData = "data p real"
            plot_cmap = 'seismic'
            u_grid.point_data[strData] = self.solution[idxF, :, 0]
        elif output == 'SPL':
            strData = "data p SPL"
            plot_cmap = 'turbo'
            u_grid.point_data[strData] = gtb.SPL(self.solution[idxF, :, 0] + self.solution[idxF, :, 1])
        else:
            Exception("output type not understood. Try 'imag', 'real' or 'SPL'.")

        u_grid.set_active_scalars(strData)
        u_plotter = pyvista.Plotter()
        u_plotter.add_mesh_isovalue(u_grid, show_edges=False, cmap=plot_cmap)
        u_plotter.add_mesh(u_grid, show_edges=False, cmap=plot_cmap, opacity='linear')
        u_plotter.view_xy()
        return u_plotter.show()

    def saveData(self, filename):
        uh_out = Function(self.V)
        with XDMFFile(self.msh.comm, '{}.xdmf'.format(filename), "w") as file:
            file.write_mesh(self.msh)
            for nf in range(len(self.freq_array)):
                uh_out.x.array.real = self.solution[nf, :, 0]
                uh_out.x.array.imag = self.solution[nf, :, 1]

                file.write_function(uh_out, self.freq_array[nf])
        return None


