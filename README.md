# electroacPy

Welcome to the electroacPy toolbox! This toolbox provides a collection of tools designed to streamline prototyping and analysis tasks in the field of electroacoustics. It consists of multiple modules, addressing specific aspects of loudspeaker system design, signal processing, and data manipulation.

## Installation

**Before "Before starting"**
These steps have been verified on Windows for Python 3.9 to 3.11. For MacOs, only 3.9 was tested. In theory, any version of Python should work as long as all dependencies are available in that specific version. Feel free to install multiple conda environment (*step 1.*) to try-out different Python version. 

**Before starting - if you don't know about installing Python**
The installation procedure is made using the CONDA package manager for Python. It can be installed either through:
- [Anaconda](https://www.anaconda.com/download/), which is a full Python development suite. It will also install Spyder (IDE), Jupyter Notebook / Lab (similar to Matlab live-script), and other Python / datascience related software, which is very practical for beginners - but a lot is a bit unnecessary. **If you chose this one, you can follow the installation guide using the newly installed Anaconda Prompt (that you'll find with the Windows / MacOs search tool).**
- [Miniconda](https://docs.anaconda.com/free/miniconda/miniconda-install/), it is "a free minimal installer for conda. It is a small bootstrap version of Anaconda that includes only conda, Python, the packages they both depend on, and a small number of other useful packages (like pip, zlib, and a few others)". **In the case you decide to use miniconda, use the Miniconda Prompt to install ElectroacPy.**
- [Miniforge](https://conda-forge.org/miniforge/), very similar to Miniconda but it is community driven (and not maintained by Anaconda) meaning you get a better architecture support for CONDA (e.g. MacOs with M1 chip). If it somewhat matters, this is the one I use. **Use Miniforge Prompt to install ElectroacPy**. Side note: for MacOs and Linux users, you'll have to install Miniforge using bash/whatever your terminal uses.

**Let's install**

1. Create a new conda environment (highly recommended but not mandatory) with the following command:

```shell
conda create -n acoustic_sim
```
   
2. Activate the environment in which electroacPy will be installed:

```shell
conda activate acoustic_sim
```

3. Set the environment to use Python 3.9 (or 3.11, see step 3+1/2). Also, if pip is not installed by default in this new environment, you can add it to the install command:

```shell
conda install python=3.9 pip
```

3+1/2. It is also possible to use Python 3.11 (and maybe Python 3.12), in that case, Tkinter (for GUI interface) comes with Python's standard library
```shell
conda install python=3.11 pip
```

4. Finally, run the following command if you don't plan contributing':

```shell
pip install /path/to/ElectroacPy
```

or, to be able to edit the toolbox:

```shell
pip install -e /path/to/ElectroacPy
```

Installing ElectroacPy in a separate environment from the "base" one will help a lot if issues are encountered when updating the toolbox. Conda environment can be selected inside the Python IDE (Spyder, PyCharm, etc.).

**Important**
If using Spyder, you'll need to install spyder-kernels in the newly created environment:
```shell
pip install spyder-kernels
```
Note that if you install spyder directly in the related conda environment, you won't have to do that.
```shell
conda install spyder
```


**Post install**
In Windows and Linux, you can actually use OpenCL to reduce computing time (and setup bempp-cl to use your computer's GPU (also works with intel CPUs)). In the corresponding Conda environment:
```shell
pip install pyopencl intel-opencl-rt
```
You might need to install OpenCL drivers, which you'll find [here](https://www.intel.com/content/www/us/en/developer/articles/technical/intel-cpu-runtime-for-opencl-applications-with-sycl-support.html). I can't remember how exactly I installed Linux opencl drivers, so you'll unfortunately have to figure that out until I look into it again.

# Modules

## electroacPy

The `electroacPy` module is a comprehensive toolkit for prototyping loudspeaker systems using both Lumped Element Method (LEM) and Boundary Element Method (BEM) approaches. It offers capabilities to design filters and crossover networks, making loudspeaker system development more efficient and effective.

## generalToolbox

The `generalToolbox` module provides a set of versatile tools that simplify various tasks related to electroacoustics and data manipulation. It offers functions to compute gains, manipulate arrays and point clouds, as well as loading UFF files (acceleration data) and directivities saved as CSV files.

## bempp-cl
From [bempp-cl website](https://bempp.com):

Bempp is an open-source computational boundary element platform to solve electrostatic, acoustic and electromagnetic problems. Features include:
- Easy-to-use Python interface.
- Support for triangular surface meshes.
- Import and export in a number of formats, including Gmsh and VTK.
- Easy formulation of acoustic and electromagnetic transmission problems.
- CPU and GPU parallelisation.
- A comprehensive operator algebra that makes it easy to formulate complex product operator formulations such as operator preconditioning.
- Coupled FEM/BEM computations via interfaces to FEniCS.

# Documentation

The PDF serves as documentation, no automatic doc was generated. 

# Contributing

If you encounter any issues or have suggestions for improvements, please feel free to contribute! You can submit issues, pull requests, or even share your usage examples in the ElectroacPy GitHub repository.

# Acknowledgments

This toolbox uses the BEMPP-CL library for Boundary Element Method computations, which is provided directly (due to some minor modifications) with electroacPy. Thanks a lot to its authors, without whom this set of modules wouldn't exists. 
