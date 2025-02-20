"""
Helpful plotting functions.

"""

import numpy as np
import electroacPy.general as gtb
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor
from matplotlib.backend_bases import MouseButton
import pyvista


#%% helper
def sumPressureArray(bemObj, radiatingSurface, radiationCoeff=None):
    p_mesh = bemObj.p_mesh
    radSurf_system = bemObj.radiatingElement
    if isinstance(radiatingSurface, str):
        if radiatingSurface == 'all':
            radiatingSurface = radSurf_system
    else:
        radiatingSurface = np.array(radiatingSurface)

    pressureCoeff = np.zeros([len(bemObj.frequency), 
                              len(bemObj.p_mesh[0, 0].coefficients)],
                             dtype='complex')
    
    if radiationCoeff is None:
        radiationCoeff = np.ones([len(bemObj.frequency), 
                                  len(radiatingSurface)], dtype=complex)
    for f in range(len(bemObj.frequency)):
        for j in range(len(radiatingSurface)):
            ind_surface = np.argwhere(radSurf_system == radiatingSurface[j])[0][0]
            pressureCoeff[f, :] += (p_mesh[f][ind_surface].coefficients * 
                                    radiationCoeff[f, j])
    return pressureCoeff


# %%2D plots
def FRF(freq, H, transformation="SPL", logx=True, legend=None, **kwargs):  #logx=True, xlim=None, ylim=None, title=None, save=False, dpi=64):
    """
    Plot Frequency Response Functions (FRFs).

    Parameters:
        freq (tuple or array): Frequency array or tuple of frequency arrays.
        H (tuple or array): FRF array or tuple of FRF arrays. If providing multiple FRF arrays,
                            each element should correspond to the FRF for the respective frequency array.
                            If a single FRF array is provided, it will be plotted against each frequency array.
        transformation (str, optional): Type of transformation to apply to the FRF data before plotting.
                                        Defaults to "SPL".
        labels (tuple or str, optional): Label(s) for the FRF(s) being plotted. If providing multiple FRFs,
                                         should be a tuple of strings or a list of strings, with each element
                                         corresponding to the label for the respective FRF(s). Defaults to None.

    Returns:
        None
    """
    from electroacPy.general.gain import SPL, dB

    # get kwargs
    labels = legend
    # check inputs
    if isinstance(freq, tuple):
        freq = freq
    else:
        freq = (freq,)

    if isinstance(H, tuple):
        H = H
    else:
        H = (H,)

    if isinstance(labels, tuple):
        labels = labels
    elif labels == None:
        pass
    else:
        labels = (labels,) * len(freq)

    if len(freq) == 1:
        if len(H) > 1:
            freq_out = freq
            for i in range(len(H)-1):
                freq_out += freq
            freq = freq_out

    # associate transformations
    if transformation in ["SPL", "spl"]:
        tr = SPL
        tr_str = "SPL [dB]"
    elif transformation == "dB":
        tr = dB
        tr_str = "gain [dB]"
    elif transformation == "abs":
        tr = np.abs
        tr_str = "magnitude [abs]"
    elif transformation == "real":
        tr = np.real
        tr_str = "Real"
    elif transformation == "imag":
        tr = np.imag
        tr_str = "Imaginary"
    elif transformation == "phase":
        tr = np.angle
        tr_str = "Phase [rad]"
    elif transformation == "unwrap_phase":
        tr = lambda x: np.unwrap(np.angle(x))
        tr_str = "Phase [rad]"
    else:
        raise Exception("transformation not understood, "
                        "available transformations: 'SPL', 'dB', 'abs', 'phase', 'real', 'imag'")

    plt.figure(figsize=kwargs.get('figsize', None))
    for i in range(len(freq)):
        if isinstance(H[i], list):
            for j in range(len(H[i])):
                if logx is True:
                    if labels == None:
                        plt.semilogx(freq[i], tr(H[i][j]))
                    else:
                        plt.semilogx(freq[i], tr(H[i][j]), label=labels[i][j])
                else:
                    if labels == None:
                        plt.plot(freq[i], tr(H[i][j]))
                    else:
                        plt.plot(freq[i], tr(H[i][j]), label=labels[i][j])
        else:
            if logx is True:
                if labels == None:
                    plt.semilogx(freq[i], tr(H[i]))
                else:
                    plt.semilogx(freq[i], tr(H[i]), label=labels[i])
            else:
                if labels == None:
                    plt.plot(freq[i], tr(H[i]))
                else:
                    plt.plot(freq[i], tr(H[i]), label=labels[i])

    plt.xlabel('Frequency [Hz]')
    # plt.title('Frequency Response Function')
    if labels != None:
        if "loc" in kwargs:
            plt.legend(loc=kwargs["loc"])
        else:
            plt.legend(loc="best")

    if 'xlim' in kwargs:
        plt.xlim(kwargs['xlim'])
    if 'ylim' in kwargs:
        plt.ylim(kwargs['ylim'])
    if 'title' in kwargs:
        plt.title(kwargs['title'])
    if 'xticks' in kwargs:
        xtick_label = []
        for i in range(len(kwargs['xticks'])):
            xtick_label.append(str(kwargs['xticks'][i]))
        plt.xticks(kwargs['xticks'], labels=xtick_label)
    if 'yticks' in kwargs:
        ytick_label = []
        for i in range(len(kwargs['yticks'])):
            ytick_label.append(str(kwargs['yticks'][i]))
        plt.yticks(kwargs['yticks'], labels=ytick_label)
    if 'ylabel' in kwargs:
        plt.ylabel(kwargs['ylabel'])
    else:
        plt.ylabel(tr_str)

    plt.grid(True, which='both', linestyle='dotted')

    if "save" in kwargs:
        if isinstance(kwargs['save'], str):
            if "dpi" in kwargs:
                plt.savefig(kwargs['save'], dpi=kwargs['dpi'])
            else:
                plt.savefig(kwargs['save'])
        else:
            print("save argument must str")
    plt.tight_layout()
    plt.show()


def directivity(freq, theta, H, transformation='SPL', logx=True, **kwargs):
    """

    :param freq:
    :param theta:
    :param H:
    :param transformation:
    :return:
    """
    from electroacPy.general.gain import dB, SPL

    # associate transformations
    if transformation == "SPL":
        tr = SPL
        tr_str = "SPL [dB]"
        tr_cmap = "turbo"
    elif transformation == "dB":
        tr = dB
        tr_str = "gain [dB]"
        tr_cmap = "turbo"
    elif transformation == "abs":
        tr = np.abs
        tr_str = "magnitude [abs]"
        tr_cmap = "turbo"
    elif transformation == "real":
        tr = np.real
        tr_str = "Real"
        tr_cmap = "seismic"
    elif transformation == "imag":
        tr = np.imag
        tr_str = "Imaginary"
        tr_cmap = "seismic"
    elif transformation == "phase":
        tr = np.angle
        tr_str = "Phase [rad]"
        tr_cmap = "magma"
    else:
        raise Exception("transformation not understood, "
                        "available transformations: 'SPL', 'dB', 'abs', 'phase', 'real', 'imag'")

    if "levels" in kwargs:
        levels = kwargs['levels']
    else:
        levels = 12

    fig, ax = plt.subplots()
    gca = ax.contourf(freq, theta, tr(H), levels, cmap=tr_cmap)
    plt.colorbar(gca, label=tr_str)
    if logx is True:
        ax.set_xscale('log')
    else:
        pass

    if 'xlim' in kwargs:
        plt.xlim(kwargs['xlim'])
    if 'ylim' in kwargs:
        plt.ylim(kwargs['ylim'])
    if 'title' in kwargs:
        plt.title(kwargs['title'])
    if 'xticks' in kwargs:
        xtick_label = []
        for i in range(len(kwargs['xticks'])):
            xtick_label.append(str(kwargs['xticks'][i]))
        plt.xticks(kwargs['xticks'], labels=xtick_label)
    if 'yticks' in kwargs:
        ytick_label = []
        for i in range(len(kwargs['yticks'])):
            ytick_label.append(str(kwargs['yticks'][i]))
        plt.yticks(kwargs['yticks'], labels=ytick_label)

    if 'ylabel' in kwargs:
        ax.set(ylabel=kwargs['ylabel'])
    else:
        ax.set(ylabel="Angle")

    # ax.set(xlabel='Frequency [Hz]', ylabel="Angle")
    ax.set(xlabel="Frequency [Hz]")
    if "save" in kwargs:
        if isinstance(kwargs['save'], str):
            if "dpi" in kwargs:
                plt.savefig(kwargs['save'], dpi=kwargs['dpi'])
            else:
                plt.savefig(kwargs['save'])
        else:
            print("save argument must str")
    plt.show()
    return None
 
def directivityViewer(theta, freq, pMic, xscale='log', fmin=20, fmax=20e3,
                         dBmin=-20, dBmax=False, title="title"):
    """
    Plot the directivity from pMic obtained with getMicPressure for a circular microphone array

    Parameters
    ----------
    theta : narray
        observation angles (in degrees).
    freq : narray
        frequency array.
    pMic : narray
        Pressur obtained with getMicPressure().
    xscale : str, optional
        scaling of the x axis (logarithmic, linear). The default is 'log'.
    fmin : float, optional
        min xlim. The default is 20.
    fmax : float, optional
        max xlim. The default is 20e3.
    dBmin : float, optional
        min ylim. The default is -20.
    dBmax : float, optional
        max ylim. The default is False.
    norm : bool, optional
        set the normalization. The default is True.

    Returns
    -------
    None.

    """
    # compute directivity
    maxAngle = np.array([np.max(abs(pMic), 1)])
    maxMatrix = np.repeat(maxAngle, len(theta), 0).T
    directivity = gtb.gain.dB(np.abs(pMic)/maxMatrix)
    dBmax = int(np.max(gtb.gain.SPL(pMic)))+3
    dBmin = dBmax-40
    
    fig = plt.figure(figsize=(11, 8))
    fig.suptitle(title)
    fig.subplots_adjust(top=0.88)
    ax1 = fig.add_subplot(221)
    gca1 = ax1.contourf(freq, theta, gtb.gain.SPL(pMic).T,
                             np.arange(dBmin, dBmax+3, 3), cmap='turbo')
    ax1.set_xscale('log')
    ax1.set(xlabel="Frequency [Hz]", ylabel="Angle [deg]", xlim=[fmin, fmax], 
            title='Raw directivity (SPL)')
    plt.colorbar(gca1)
    
    # Total directivity NORMALIZED
    ax2 = fig.add_subplot(223)
    gca2 = ax2.contourf(freq, theta, directivity.T, np.arange(-21, 3, 3), cmap='turbo')
    ax2.set_xscale('log')
    ax2.set(xlabel="Frequency [Hz]", ylabel="Angle [deg]", xlim=[fmin, fmax], 
            title='Absolute directivity')
    plt.colorbar(gca2)
    
    # pressure response
    obsAngle = 0
    ind_angle, value_angle = gtb.findInArray(theta, np.deg2rad(obsAngle))
    ax3 = fig.add_subplot(222)
    ax3.semilogx(freq, gtb.gain.SPL(pMic[:, ind_angle]), label=int(np.rad2deg(value_angle)))
    ax3.set(xlabel='Frequency [Hz]', ylabel='SPL [dB]', ylim=[dBmin, dBmax], 
            title='Pressure response')
    ax3.grid(which='both')
    ax3.legend(loc='best')
    
    # polar
    freqObs = 1000
    ind_freq, value_freq = gtb.findInArray(freq, freqObs)
    ax4 = fig.add_subplot(224, polar=True)
    ax4.plot(np.deg2rad(theta), gtb.gain.SPL(pMic[ind_freq, :]), label=int(value_freq))
    ax4.set(ylim=[dBmin, dBmax], title='Polar directivity')
    ax4.legend(loc='best')

    #defining the cursor
    cursor = Cursor(ax1, horizOn = True, vertOn=True, color='black', linewidth=1.2, 
                    useblit=True)
    
    def onclick(event):
        if event.button is MouseButton.LEFT:
            x = event.xdata
            y = event.ydata
             
            # update subplots -- FRF
            obsAngle = y
            ind_angle, value_angle = gtb.findInArray(theta, obsAngle)
            ax3.semilogx(freq, gtb.gain.SPL(pMic[:, ind_angle]), label=int(value_angle))
            ax3.legend(loc='best')
            
            # update subplots -- polar
            freqObs = x
            ind_freq, value_freq = gtb.findInArray(freq, freqObs)
            ax4.plot(np.deg2rad(theta), 
                     gtb.gain.SPL(pMic[ind_freq, :]), label=int(value_freq))
            ax4.legend(loc='best')
            
            fig.canvas.draw() #redraw the figure
        
        if event.button is MouseButton.RIGHT:
            if len(ax3.lines) > 1:
                ax3.lines[-1].remove()
                ax4.lines[-1].remove()
                ax3.legend(loc='best')
                ax4.legend(loc='best')
                fig.canvas.draw()
        
    fig.canvas.mpl_connect('button_press_event', onclick)
    plt.tight_layout()
    plt.show()
    return cursor
 
    
#%% 3D plots
def bempp_grid(bemOBJ, eval_grid, pmic, radiationCoeff, radiatingElement,
               eval_name, transformation="SPL", export_grid=False):
    """
    Plot grid functions from bempp. "grids" argument takes a list of one or
    multiple grids (evaluation grids), and concatenates everything based on 
    which radiating elements are selected.

    Parameters
    ----------
    bemOBJ : bem object
        bempp simulation,
    eval_grid : list
        eval_grid to plot.
    p_mic : ndarray
        pressure received at microphones (one per eval_grid[i])
    radiationCoeff : ndarray
        coefficients to apply on surfaces
    radiatingElement : list
        radiating elements to sum on system mesh surface.
    eval_name: list of str
        names of eval_grid components
    transformation : str
        scale plot to "SPL" (or "spl"), "real", "imag" or "phase".

    Returns
    -------
    GMSH or Paraview plot.

    """
    import gmsh
    import bempp.api
    import tempfile
    import subprocess

    # get transform
    if transformation in ["SPL", "spl"]:
        T = gtb.gain.SPL
    elif transformation == "real":
        T = np.real
    elif transformation == "imag":
        T = np.imag
    elif transformation == "phase":
        T = np.angle
    
    # get some data from system mesh
    grid         = bemOBJ.grid_init
    Nvert_grid   = grid.vertices.shape[1]
    frequency    = bemOBJ.frequency
    
    # get pressure_grid
    pressure_grid = sumPressureArray(bemOBJ, radiatingElement, radiationCoeff)
    pressure_grid = pressure_grid[:, :Nvert_grid]
          
    pressure = {"system": pressure_grid}
    pressure_id = {"system": 0}
    for i in range(len(eval_name)):
        pressure[eval_name[i]] = pmic[i]
        pressure_id[eval_name[i]] = i+1
        
    # set all grid and plotting grid to a single list (for union!)
    grids = []
    grids.append(bempp.api.Grid(grid.vertices, grid.elements))
    coeff_limits = np.zeros(len(eval_grid)+2, dtype=int)
    coeff_limits[0] = 1
    coeff_limits[1] = Nvert_grid+1
    for i in range(len(eval_name)):
        grid_tmp = bempp.api.Grid(eval_grid[i].vertices, eval_grid[i].elements)
        grids.append(grid_tmp)
        coeff_limits[i+2] = coeff_limits[i+1]+eval_grid[i].vertices.shape[1]

    grid_union = bempp.api.grid.union(grids)
    
    with tempfile.NamedTemporaryFile(delete=False, suffix=".msh") as tmp_file:
        bempp.api.export(tmp_file.name, grid_union, write_binary=False)

    
    # Initialize GMSH
    gmsh.initialize()
    gmsh.option.setNumber("PostProcessing.Format", 5)
    gmsh.option.setNumber("PostProcessing.Link", 1)
    
    # Load total grid and assign view, step, value, etc.
    gmsh.open(tmp_file.name)
    view_tag = gmsh.view.add("evaluation")
    model = gmsh.model.getCurrent()
    for i, surface in enumerate(gmsh.model.getEntities(2)):  # 2 = surfaces
        gmsh.model.addPhysicalGroup(2, [surface[1]], tag=i)
    
    for ps in pressure: 
        ps_id = pressure_id[ps]
        N1 = coeff_limits[ps_id] 
        N2 = coeff_limits[ps_id+1] 
        for j, f in enumerate(frequency):
            data = T(pressure[ps][j, :])
            node_tags = np.arange(N1, N2)
            gmsh.view.addModelData(view_tag, 
                                   j, 
                                   model,
                                   "NodeData", 
                                   node_tags,
                                   data[:, None], 
                                   time=f)
    
    with tempfile.NamedTemporaryFile(delete=False, suffix=".msh") as tmp_export:
        gmsh.write(tmp_export.name)
        gmsh.view.write(view_tag, tmp_export.name)
        
    if export_grid is not False:
        gmsh.write(export_grid)
        gmsh.view.write(view_tag, export_grid)
    
    gmsh.finalize()
    
    # Open the combined file in GMSH viewer
    subprocess.Popen(["gmsh", 
                      "-setnumber", "Mesh.SurfaceEdges", "0", 
                      tmp_export.name], 
                     stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
    return None

def bempp_grid_mesh(bemOBJ, radiationCoeff, radiatingElement, 
                  transformation="SPL", export_grid=False):
    """
    Plot grid functions from bempp. "grids" argument takes a list of one or
    multiple grids (evaluation grids), and concatenates everything based on 
    which radiating elements are selected.

    Parameters
    ----------
    bemOBJ : bem object
        bempp simulation,
    radiationCoeff : ndarray
        coefficients to apply on surfaces
    radiatingElement : list
        radiating elements to sum on system mesh surface.
    transformation : str
        scale plot to "SPL" (or "spl"), "real", "imag" or "phase".

    Returns
    -------
    GMSH or Paraview plot.

    """
    import gmsh
    import tempfile
    import subprocess

    # get transform
    if transformation in ["SPL", "spl"]:
        T = gtb.gain.SPL
    elif transformation == "real":
        T = np.real
    elif transformation == "imag":
        T = np.imag
    elif transformation == "phase":
        T = np.angle
    
    # get some data from system mesh
    grid         = bemOBJ.grid_init
    Nvert_grid   = grid.vertices.shape[1]
    frequency    = bemOBJ.frequency
    
    # get pressure_grid
    pressure_grid = sumPressureArray(bemOBJ, radiatingElement, radiationCoeff)
    pressure_grid = pressure_grid[:, :Nvert_grid]          
    pressure = {"system": pressure_grid}
        
    # Initialize GMSH
    gmsh.initialize()
    gmsh.option.setNumber("PostProcessing.Format", 5)
    gmsh.option.setNumber("PostProcessing.Link", 1)
    
    # Load total grid and assign view, step, value, etc.
    gmsh.open(bemOBJ.meshPath)
    view_tag = gmsh.view.add("system pressure")
    model = gmsh.model.getCurrent()
    for i, surface in enumerate(gmsh.model.getEntities(2)):
        gmsh.model.addPhysicalGroup(2, [surface[1]], tag=i)
    
    node_tags = np.arange(1, len(bemOBJ.p_mesh[0, 0].coefficients)+1)
    for ps in pressure: 
        for j, f in enumerate(frequency):
            data = T(pressure[ps][j, :])
            gmsh.view.addModelData(view_tag, 
                                   j, 
                                   model,
                                   "NodeData", 
                                   node_tags,
                                   data[:, None], 
                                   time=f)
    
    with tempfile.NamedTemporaryFile(delete=False, suffix=".msh") as tmp_export:
        gmsh.write(tmp_export.name)
        gmsh.view.write(view_tag, tmp_export.name)
        
    if export_grid is not False:
        gmsh.write(export_grid)
        gmsh.view.write(view_tag, export_grid)
    
    gmsh.finalize()
    
    # Open the combined file in GMSH viewer
    subprocess.Popen(["gmsh", 
                      "-setnumber", "Mesh.SurfaceEdges", "0", 
                      tmp_export.name], 
                     stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
    return None

def pressureField_3D(bemOBJ, xMic, L, W, pMicData, radiatingElement):
    """
    Plot pressure field using PyVista backend.
    :param bemOBJ: object,
        bem object used in the current observation.
    :param L: list,
        list of length of each planar surface defined using "planar" observations.
    :param W: list,
        list of width of each planar surface defined using "planar" observations.
    :param pMicData: ndarray,
        pressure field data.
    :param xMic: ndarray,
        microphones position.
    :param radiatingElement: str, list,
        Index of surfaces currently radiating.
    :return:
    """
    # initialize inputs
    freq_array = bemOBJ.frequency
    meshPlot = []
    
    # put inputs into list (facilitate manipulation when switching between
    # multiple/single plot)
    nObs = len(xMic)
    
    # initialize plotter
    if bemOBJ.domain == "interior":
        pl = pyvista.Plotter(lighting='three lights', shape=(1, 2))
    elif bemOBJ.domain =="exterior":
        pl = pyvista.Plotter(lighting='three lights')
    pl.background_color = 'white'

    for i in range(nObs):
        meshPlot.append(pyvista.StructuredGrid())
        
        # initialisation
        meshPlot[i].points = xMic[i]
        meshPlot[i].points = xMic[i]
        meshPlot[i].dimensions = [len(L[i]), len(W[i]), 1]


    engine = pyvista_spkPlotter(pl, bemOBJ, L, W,
                                pMicData, xMic, radiatingElement, nObs)

    def plotPressurePoint(point):
        position = np.zeros([nObs, 3])
        for i in range(len(xMic)):
            _, position[i, :] = gtb.geometry.findClosestPoint(xMic[i], point[0], point[1], point[2])
        idxp = gtb.geometry.findClosestPoint(position, point[0], point[1], point[2])[0]
        idxm, p = gtb.geometry.findClosestPoint(xMic[idxp], point[0], point[1], point[2])

        pMic_tmp = pMicData[idxp].T
        pMicToPlot = pMic_tmp[idxm, :]
        maxP = np.max(gtb.gain.SPL(pMicToPlot))

        fig, ax = plt.subplots()
        ax.semilogx(bemOBJ.frequency, gtb.gain.SPL(pMicToPlot), label='{}'.format(p))
        ax.grid(which='both')
        ax.set(xlabel='Frequency [Hz]', ylabel='SPL [dB]', ylim=[maxP-40, 
                                                                 maxP+6])
        ax.legend(loc='best')
        plt.show()
        return None


    pl.add_slider_widget(
        callback=lambda value: engine('cfreq', value),
        rng=[np.log10(freq_array[0]), np.log10(freq_array[-1])],
        value=np.log10(250),
        title="Frequency [log]",
        style="modern"
    )

    pl.add_checkbox_button_widget(
        callback=lambda value: engine('real', value),
        value=False,
        position=(5, 100),
    )

    pl.add_checkbox_button_widget(
        callback=lambda value:engine('showMesh', value),
        value=False,
        position=(5, 0),
    )

    pl.add_checkbox_button_widget(
        callback=lambda value:engine('contourP', value),
        value=True,
        position=(5, 50),
    )

    pl.enable_surface_point_picking(
        callback=plotPressurePoint,
        show_message=False,
    )

    obj = pl.show()
    return obj


def boundingBox(bemOBJ, nx, ny, nz, pMicData, xMic, radiatingElement):
    """
    Plot a bounding box observation in a PyVista Plotter
    :param bemOBJ: object,
        bem object used in the current observation.
    :param pMicData: ndarray,
        pressure field data.
    :param xMic: ndarray,
        microphones position.
    :param radiatingElement: str, list,
        Index of surfaces currently radiating.
    :return:
    """
    freq_array = bemOBJ.frequency
    meshPath = bemOBJ.meshPath
    meshPlot = []

    # init plotter
    pl = pyvista.Plotter(lighting='three lights')
    pl.background_color = 'white'
    sizeFactor = bemOBJ.sizeFactor
    vertices = bemOBJ.vertices

    engine = pyvista_boundingBoxPlotter(pl, bemOBJ, nx, ny, nz, pMicData, xMic, radiatingElement)

    pl.add_slider_widget(
        callback=lambda value:engine('cfreq', value),
        rng=[np.log10(freq_array[0]), np.log10(freq_array[-1])],
        value=np.log10(250),
        title="Frequency, log()",
        style="modern"
    )

    pl.add_checkbox_button_widget(
        callback=lambda value:engine('real', value),
        value=False,
        position=(5, 100),
    )

    pl.add_checkbox_button_widget(
        callback=lambda value:engine('showMesh', value),
        value=False,
        position=(5, 0),
    )

    pl.add_checkbox_button_widget(
        callback=lambda value:engine('contourP', value),
        value=True,
        position=(5, 50),
    )

    pl.add_plane_widget(
        callback=lambda normal, origin: engine('normOrigin', (normal, origin)),
        normal=(1, 0, 0),
        origin=(0, 0, 0),
        normal_rotation=True,
    )

    obj = pl.show()
    return obj

def sphericalRadiation(xMic, pMic, freq_array):
    # init plotter
    pl = pyvista.Plotter(lighting='three lights')
    pl.background_color = 'white'

    engine = pyvista_sphericalPlotter(pl, pMic, xMic, freq_array)

    pl.add_slider_widget(
        callback=lambda value:engine('cfreq', value),
        rng=[np.log10(freq_array[0]), np.log10(freq_array[-1])],
        value=np.log10(250),
        title="Frequency, log()",
        style="modern"
    )
    obj = pl.show()
    return obj

## ================================================
# %% PyVista plotter classes (for 3D visualization)
class pyvista_spkPlotter:
    # I think this class and associated methods are more complex than 
    # they should be. For now it works as it is.
    def __init__(self, mesh, bemOBJ, L, W, pMicData, 
                 xMic, radiatingElement, nObs):
        self.output = mesh
        self.bemOBJ = bemOBJ
        self.freq_array = bemOBJ.frequency
        self.vertices = bemOBJ.vertices
        self.sizeFactor = bemOBJ.sizeFactor
        self.L = L
        self.W = W
        self.pMicData = pMicData
        self.xMic = xMic
        self.radiatingElement = radiatingElement
        self.nObs = len(pMicData)

        self.speaker = pyvista.read(bemOBJ.meshPath)
        meshPlot = []
        # pmicdata
        for i in range(nObs):
            meshPlot.append(pyvista.StructuredGrid())
            # initialisation
            meshPlot[i].points = xMic[i]
            meshPlot[i].dimensions = [len(L[i]), len(W[i]), 1]
        self.meshPlot = meshPlot

        # default parameters
        self.kwargs = {
            'real': False,
            'cfreq': '250',
            'showMesh': False,
            'contourP': True
        }

    def __call__(self, param, value):
        self.kwargs[param] = value
        self.update()

    def update(self):
        real = self.kwargs['real']
        cfreq = self.kwargs['cfreq']
        showMesh = self.kwargs['showMesh']
        contourP = self.kwargs['contourP']

        # min / max bounds - title
        maxPressure = []
        minPressure = []
        if real is True:
            for i in range(self.nObs):
                maxPressure.append(np.max(np.real(self.pMicData[i])))
                minPressure.append(np.min(np.real(self.pMicData[i])))
            minPressure = -1 #np.min(minPressure)
            maxPressure = 1  #np.max(maxPressure)
            pressureCoeff_system = sumPressureArray(self.bemOBJ, 
                                                    self.radiatingElement)
            pressureCoeff_system = np.real(pressureCoeff_system)/10
            title_arg = 'Real(Pressure), normalized - '
            cmap_d = 'seismic'
        else:
            for i in range(self.nObs):
                maxPressure.append(np.max(gtb.gain.SPL(self.pMicData[i])))

            maxPressure = np.max(maxPressure) + 5
            minPressure = maxPressure - 75

            # to normalize pressureCoeff_system
            # (pressure on baffle will usually be super high compared to radiated field)
            maxP_abs = 10**(maxPressure/20) * 2e-5

            pressureCoeff_system = sumPressureArray(self.bemOBJ, np.array(self.radiatingElement))
            pressureCoeff_system = gtb.normMinMax(np.abs(pressureCoeff_system), 0, maxP_abs)
            pressureCoeff_system = gtb.gain.SPL(pressureCoeff_system)
            title_arg = 'SPL, dB - '
            cmap_d = 'turbo'

        if contourP is True:
            ncolor=21
        else:
            ncolor=256

        freq_idx, frequency = gtb.findInArray(self.freq_array, 10 ** (cfreq))
        sargs = dict(
            title_font_size=21,
            label_font_size=18,
            shadow=False,
            n_labels=5,
            italic=False,
            font_family="arial",
            color='k',
            title=title_arg + str(round(frequency, 1)) + " Hz",
        )
        
        # system mesh
        try:
            if self.bemOBJ.domain == "exterior":
                self.output.remove_actor('spk_mesh')
                self.output.add_mesh(self.speaker, show_edges=showMesh,
                            scalars=pressureCoeff_system[freq_idx, :self.bemOBJ.spaceP.grid_dof_count // self.sizeFactor],
                            cmap=cmap_d,
                            n_colors=ncolor,
                            show_scalar_bar=False, clim=[minPressure, maxPressure], name='spk_mesh')

            elif self.bemOBJ.domain == "interior":
                self.output.remove_actor('spk_mesh')
                self.output.subplot(0, 0)
                self.output.add_mesh(self.speaker, show_edges=showMesh,
                            scalars=pressureCoeff_system[freq_idx, :self.bemOBJ.spaceP.grid_dof_count // self.sizeFactor],
                            cmap=cmap_d,
                            n_colors=ncolor,
                            show_scalar_bar=False, clim=[minPressure, maxPressure], 
                            name='spk_mesh')
        except:
            pressure_system = pressureCoeff_system[freq_idx, :self.bemOBJ.spaceP.grid_dof_count // self.sizeFactor]
            toAdd = self.vertices - self.bemOBJ.spaceP.grid_dof_count
            scalars_to_plot = np.concatenate((pressure_system, np.zeros(toAdd)))
            if self.bemOBJ.domain == "exterior":
                self.output.remove_actor('spk_mesh')
                self.output.add_mesh(self.speaker, show_edges=showMesh,
                                     scalars=scalars_to_plot,
                                     cmap=cmap_d,
                                     n_colors=ncolor,
                                     show_scalar_bar=False, clim=[minPressure, maxPressure],
                                     name='spk_mesh')
            elif self.bemOBJ.domain == "interior":
                self.output.remove_actor('spk_mesh')
                self.output.subplot(0, 0)
                self.output.add_mesh(self.speaker, show_edges=showMesh,
                            scalars=pressureCoeff_system[freq_idx, :self.bemOBJ.spaceP.grid_dof_count // self.sizeFactor],
                            cmap=cmap_d,
                            n_colors=ncolor,
                            show_scalar_bar=False, 
                            clim=[minPressure, maxPressure], 
                            name='spk_mesh')
        
        show_sc_b = []
        for i in range(self.nObs):
            if i == 0:
                show_sc_b.append(True)
            else:
                show_sc_b.append(False)

            data = self.pMicData[i].T
            pMicRes = np.reshape(data[:, int(freq_idx)],
                                 [len(self.L[i]), len(self.W[i])])

            # def real or SPL pressure plot
            if real is True:
                # pMicToPlot = gtb.normMinMax(np.real(pMicRes), minPressure, maxPressure)
                pMicToPlot = np.real(pMicRes)
            else:
                pMicToPlot = gtb.gain.SPL(pMicRes)

            # observation mesh
            if self.bemOBJ.domain == "exterior":
                self.output.remove_actor('SPL_{}'.format(i))
                self.output.add_mesh(self.meshPlot[i], show_edges=showMesh, scalars=pMicToPlot, cmap=cmap_d,
                            scalar_bar_args=sargs, clim=[minPressure, maxPressure], n_colors=ncolor,
                            name='SPL_{}'.format(i), show_scalar_bar=show_sc_b[i])
            elif self.bemOBJ.domain == "interior":
                self.output.subplot(0, 1)
                self.output.remove_actor('SPL_{}'.format(i))
                self.output.add_mesh(self.meshPlot[i], show_edges=showMesh, scalars=pMicToPlot, cmap=cmap_d,
                            scalar_bar_args=sargs, clim=[minPressure, maxPressure], n_colors=ncolor,
                            name='SPL_{}'.format(i), show_scalar_bar=show_sc_b[i])
                self.output.add_mesh(self.speaker, show_scalar_bar=False, color='grey',
                                     opacity=0.15, name='boundary_{}'.format(i))
                self.output.link_views()
        
        if self.bemOBJ.domain == "exterior": # if used for interior plot, crashes
            self.output.add_axes(color='k') # this mf crashes
        
        _ = self.output.show_grid(color='k') # this one does not
        return


class pyvista_boundingBoxPlotter:
    def __init__(self, mesh, bemOBJ, nx, ny, nz, pMicData, xMic, radiatingElement):
        self.output = mesh
        self.bemOBJ = bemOBJ
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.freq_array = bemOBJ.frequency
        self.vertices = bemOBJ.vertices
        self.sizeFactor = bemOBJ.sizeFactor
        self.pMicData = pMicData
        self.xMic = xMic
        self.radiatingElement = radiatingElement

        self.speaker = pyvista.read(bemOBJ.meshPath) # will display the speaker's mesh

        # pmicdata
        meshPlot = pyvista.StructuredGrid()
        meshPlot.points = xMic
        meshPlot.dimensions = [nx, ny, nz]
        self.meshPlot = meshPlot   # will display bounding box mesh (unstructured grid?)

        # default parameters
        self.kwargs = {
            'real': False,
            'cfreq': '250',
            'showMesh': False,
            'contourP': True,
            'normal': (1, 0, 0),
            'origin': (0, 0, 0)
        }

    def __call__(self, param, value):
        self.kwargs[param] = value
        self.update()

    def update(self):
        bemOBJ = self.bemOBJ
        real = self.kwargs['real']
        cfreq = self.kwargs['cfreq']
        showMesh = self.kwargs['showMesh']
        contourP = self.kwargs['contourP']
        # normOrigin = self.kwargs['normOrigin']
        normal = self.kwargs['normal']
        origin = self.kwargs['origin']

    
        # normal = normOrigin[0]
        # origin = normOrigin[1]

        # min / max bounds - title
        if real is True:
            minPressure = np.min(np.real(self.pMicData))
            maxPressure = np.max(np.real(self.pMicData))
            pressureCoeff_system = sumPressureArray(self.bemOBJ, self.radiatingElement)
            pressureCoeff_system = gtb.normMinMax(np.real(pressureCoeff_system), minPressure, maxPressure)
            title_arg = 'Real(Pressure), Pa - '
            cmap_d = 'seismic'
        else:
            maxPressure = np.max(gtb.gain.SPL(self.pMicData)) + 5
            minPressure = maxPressure - 75

            # to normalize pressureCoeff_system
            # (pressure on baffle will usually be super high compared to radiated field)
            maxP_abs = 10 ** (maxPressure / 20) * 2e-5

            pressureCoeff_system = sumPressureArray(self.bemOBJ, self.radiatingElement)
            pressureCoeff_system = gtb.normMinMax(np.abs(pressureCoeff_system), 0, maxP_abs)
            pressureCoeff_system = gtb.gain.SPL(pressureCoeff_system)
            title_arg = 'SPL, dB - '
            cmap_d = 'turbo'

        if contourP is True:
            ncolor = 21
        else:
            ncolor = 256

        freq_idx, frequency = gtb.findInArray(self.freq_array, 10 ** (cfreq))
        sargs = dict(
            title_font_size=21,
            label_font_size=18,
            shadow=False,
            n_labels=5,
            italic=False,
            font_family="arial",
            color='k',
            title=title_arg + str(round(frequency, 1)) + " Hz",
        )


        data = self.pMicData.T
        pMicRes = data[:, int(freq_idx)]

        # def real or SPL pressure plot
        if real is True:
            pMicToPlot = gtb.normMinMax(np.real(pMicRes), minPressure, maxPressure)
        else:
            pMicToPlot = gtb.gain.SPL(pMicRes)

        # system mesh
        _ = self.output.remove_actor('spk_mesh')
        self.output.add_mesh(self.speaker, show_edges=showMesh,
                             scalars=pressureCoeff_system[freq_idx, :self.vertices // self.sizeFactor],
                             cmap=cmap_d, n_colors=ncolor, show_scalar_bar=False,
                             clim=[minPressure, maxPressure], name='spk_mesh')

        # observation mesh
        meshPlot = self.meshPlot
        meshPlot.point_data['scalars'] = pMicToPlot
        slc = meshPlot.slice(normal=normal, origin=origin)
        # self.output.add_mesh(meshPlot.outline(), color='k')
        self.output.add_mesh(slc, show_edges=showMesh, show_scalar_bar=True,
                                   cmap=cmap_d,  clim=[minPressure, maxPressure], n_colors=ncolor,
                                   name='SPL_slice', scalar_bar_args=sargs)
        _ = self.output.show_grid(color='k')
        self.output.add_axes(color='k')
        return


class pyvista_sphericalPlotter:
    def __init__(self, mesh, pMic, xMic, frequency):
        self.output = mesh
        self.pMic = pMic
        self.xMic = xMic
        self.frequency = frequency

        ELEVATION = np.arccos(xMic[:, 2] / np.sqrt(xMic[:, 0] ** 2 + xMic[:, 1] ** 2 + xMic[:, 2] ** 2))
        AZIMUTH = np.zeros(len(xMic))
        self.R = np.zeros([len(frequency), len(xMic)])
        for i in range(len(xMic)):
            if xMic[i, 1] > 0:
                AZIMUTH[i] = np.arccos(xMic[i, 0] / np.sqrt(xMic[i, 0] ** 2 + xMic[i, 1] ** 2))
            elif xMic[i, 1] < 0:
                AZIMUTH[i] = - np.arccos(xMic[i, 0] / np.sqrt(xMic[i, 0] ** 2 + xMic[i, 1] ** 2))

        self.elevation = ELEVATION
        self.azimuth = AZIMUTH

        self.kwargs = {
            'cfreq': '250',
        }

    def __call__(self, param, value):
        self.kwargs[param] = value
        self.update()

    def update(self):
        cfreq = self.kwargs['cfreq']
        # print(cfreq)
        R = self.R
        elevation = self.elevation
        azimuth = self.azimuth

        freq_idx, frequency = gtb.findInArray(self.frequency, 10 ** (cfreq))

        sargs = dict(
            title_font_size=21,
            label_font_size=18,
            shadow=False,
            n_labels=5,
            italic=False,
            font_family="arial",
            color='k',
            title="SPL, dB - " + str(round(frequency, 1)) + " Hz",
        )

        cloud = pyvista.PolyData(self.xMic)
        volume = cloud.delaunay_3d(alpha=2.)
        shell = volume.extract_geometry()
        shell.points /= np.sqrt(self.xMic[0, 0]**2 + self.xMic[0, 1]**2 + self.xMic[0, 2]**2)
        shell.points *= np.tile(gtb.gain.SPL(self.pMic[freq_idx, :]), (3, 1)).T

        _ = self.output.remove_actor('balloon')
        self.output.add_mesh(shell, name='balloon', cmap="turbo",
                             scalars=gtb.gain.SPL(self.pMic[freq_idx, :]),
                             clim=[np.max(gtb.gain.SPL(self.pMic))-30, np.max(gtb.gain.SPL(self.pMic)) + 6],
                             n_colors=20, scalar_bar_args=sargs)
        _ = self.output.show_grid(color='k')
        self.output.add_axes(color='k')
        return None


