# electroacPy

Welcome to the electroacPy toolbox! This toolbox provides a collection of tools designed to streamline prototyping and analysis tasks in the field of electroacoustics. It consists of multiple modules, addressing specific aspects of loudspeaker system design, signal processing, and data manipulation.

## Installation
**Before Starting**
- The following steps have been verified on **Windows** for Python versions 3.9 to 3.11. For **macOS**, only version 3.9 has been tested. However, in theory, any version of Python should work as long as all dependencies are available.
- You may want to try out different Python versions by creating multiple Conda environments (see **Step 1**).

### Setting Up Python with Conda
The recommended installation method uses the **Conda** package manager for Python. You can install Conda through one of the following options:
1. [Anaconda](https://www.anaconda.com/download/): A full Python development suite that includes Spyder (IDE), Jupyter Notebook/Lab, and other tools.
    - **Windows**: Use the Anaconda Prompt to follow the installation steps.
    - **macOS/Linux**: Use your terminal (bash/zsh).
2. [Miniconda](https://docs.anaconda.com/free/miniconda/miniconda-install/): A minimal version of Anaconda, including only the necessary packages for managing environments.
    - **Windows**: Use the Miniconda Prompt for installation.
    - **macOS/Linux**: Use your terminal (bash/zsh).
3. [Miniforge](https://conda-forge.org/miniforge/): A lightweight version similar to Miniconda, but community-driven, providing better architecture support (e.g., M1/M2 chips on macOS).
    - **Windows**: Use the Miniforge Prompt for installation.
    - **macOS/Linux**: Use your terminal (bash/zsh).

### Installation Steps
1. **Create a new Conda environment** (recommended but optional):
```bash
conda create -n acoustic_sim
```

2. **Activate the environment**:
```bash
conda activate acoustic_sim
```

3. **Install Python 3.9 and pip** (you can adjust the Python version if needed):
```bash
conda install python=3.9 pip
```
Optionally, you can use Python 3.11 (or newer).
```bash
conda install python=3.11
```
4. **Install electroacPy**:
- For standard installation:
    ```bash
    pip install /path/to/electroacPy
    ```
- For development (editable) installation:
    ```bash
    pip install -e /path/to/electroacPy
    ```

### Notes
- **Using a separate environment**:  Installing ElectroacPy in its own Conda environment is recommended. This helps prevent conflicts during updates and allows easier management of dependencies.
- **Selecting environments**: In Python IDEs like Spyder or PyCharm, you can choose the specific Conda environment where ElectroacPy is installed.

## Additional Steps for Spyder Users
If you plan to use **Spyder**:
- You'll need to install `spyder-kernels` in the newly created environment:
```shell
pip install spyder-kernels
```
- Alternatively, you can install **Spyder** directly in the environment to avoid needing `spyder-kernels`:
```shell
conda install spyder
```

## Post install
In Windows and Linux, you can actually use the OpenCL backend to reduce computing time. In the corresponding Conda environment:
```shell
pip install pyopencl intel-opencl-rt
```
You'll also need to install OpenCL drivers, which you'll find [here](https://www.intel.com/content/www/us/en/developer/articles/technical/intel-cpu-runtime-for-opencl-applications-with-sycl-support.html) for windows users. Linux users can follow the **OpenCL** section from [bempp-cl installation guide](https://bempp.com/installation.html).

# Modules
This toolbox is divided into multiple sub-modules.

## electroacPy

The `electroacPy` module is a toolkit for prototyping loudspeaker systems using Lumped Element Method (LEM). It is also a set of wrappers for bempp-cl to solve acoustic radiation problems using Boundary Elements. It offers capabilities to design filters and crossover networks.

## generalToolbox

The `generalToolbox` module provides a set of tools that simplify various tasks related to electroacoustics and data manipulation. It offers functions to compute gains, manipulate arrays and point clouds, plot frequency-response functions, as well as loading UFF files (acceleration data) and directivities saved as CSV files.

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

Waiting for documentation. For now, please refer to the examples.

# Contributing

If you encounter any issues or have suggestions for improvements, please feel free to fork and contribute. You can submit issues, pull requests, or even share your usage examples.

# Acknowledgments

This toolbox uses the BEMPP-CL library for Boundary Element Method computations, which is provided directly (due to some minor modifications) in this repo. Thanks a lot to its authors. 
