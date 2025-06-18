# electroacPy

This module provides a collection of tools designed to streamline prototyping and analysis tasks in the field of electroacoustics. It consists of multiple modules, addressing specific aspects of loudspeaker-system design.

The main GitHub repository is located [here](https://github.com/TomMunoz3772/electroacPy).
Examples are available in [another repository](https://github.com/TomMunoz3772/electroacPy_examples).
Documentation of the main repository is available on [readthedocs](https://electroacpy.readthedocs.io).

## Installation
### Before Starting

- The following steps have been verified on **Windows** and **Linux** for Python versions 3.9 to 3.12. For **macOS**, only version 3.9 has been tested. However, in theory, any version of Python should work as long as all dependencies are available.
- You may want to try out different Python versions by creating multiple Conda environments.


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

#### Creating a conda environment
If you decide to manage your Python environments with conda, you can setup a specific environment for electroacPy as follow:

1. **Create the environment**
```shell
conda create -n acoustic_sim
```

2. **Activate the environment**
```shell
conda activate acoustic_sim
```

3. **Install a version of Python compatible with electroacPy, as well as the *pip* package manager**
```shell
conda install python=3.12 pip
```

These 3 steps are all you need to do before installing electroacPy.

### Recommended installation: PyPI
The easiest way to install electroacpy is to use the [PyPI repository](https://pypi.org/project/electroacPy/). In your Python/Conda environment:
```shell
pip install electroacPy
```

### Install from source
If however, you wish to install from source, you'll need first to download or clone the main repository --- [here's a link to help you understand the process](https://docs.github.com/en/get-started/start-your-journey/downloading-files-from-github). Once you have downloaded what you need, you can use the following command (remember to be in your Python/Conda environment):

```shell
pip install /path/to/electroacPy
```

For development installation:
```shell
pip install -e /path/to/electroacPy
```

You'll need to replace `/path/to/electroacPy` by the actual location of  the toolbox on your computer --- pointing to the folder containing the "pyproject.toml" file. For example, if you use Windows, the path can look like this: `C:\Users\yourUsername\Documents\GitHub\electroacPy`.

### Notes
**Using a separate environment**:  Installing ElectroacPy in its own Conda / Python environment is recommended. This helps prevent conflicts during updates and allows easier management of dependencies.

**Selecting environments**: In Python IDEs like Spyder or PyCharm, you can choose the specific Conda environment where ElectroacPy is installed.

### Additional Steps for Spyder Users
If you plan to use **Spyder**:

You'll need to install `spyder-kernels` in the newly created environment:
```shell
pip install spyder-kernels
```

Alternatively, you can install **Spyder** directly in the environment to avoid needing `spyder-kernels`:
```shell
conda install spyder
```

## OpenCL
In Windows and Linux, you can actually use the OpenCL backend to reduce computing time. In the corresponding Conda environment:
```shell
pip install pyopencl
```

You'll also need to install OpenCL drivers, which you'll find [here](https://www.intel.com/content/www/us/en/developer/articles/technical/intel-cpu-runtime-for-opencl-applications-with-sycl-support.html) for intel users. For more information, you can follow the **OpenCL** section from [bempp-cl installation guide](https://bempp.com/installation.html).

---

# Modules
## electroacPy

The `ElectroacPy` module is a toolkit for prototyping loudspeaker systems using Lumped Element Method (LEM). It is also a set of wrappers for bempp-cl to solve acoustic radiation problems using Boundary Elements. It offers capabilities to design filters and crossover networks.

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

ElectroacPy's handbook can be found [here](https://electroacpy.readthedocs.io). Full documentation is still in progress.

# Contributing

If you encounter any issues or have suggestions for improvements, please feel free to fork and contribute. You can submit issues, pull requests, or even share your usage examples.

# Acknowledgments

This toolbox uses the [bempp-cl](https://bempp.com) library for Boundary Element Method computations, which is provided directly in this repo. Many thanks to its authors and contributors, without whom this toolbox wouldn't have been available in its current form.

# Licence
This toolbox, ElectroacPy, is licensed under the GNU General Public License, Version 3 (GPLv3).

Certain parts of the source code use bempp-cl, which is licensed under the MIT License. The MIT-licensed portions remain under their original license.

