# Introduction
The main idea behind this project is to provide an easy-to-use toolbox to solve simple acoustic studies using Python --- with focus on loudspeaker radiation. It uses a mix of lumped-element-modeling (LEM) for drivers and boundary-element method (BEM) for interior and exterior acoustical radiation. This handbook provides examples on how to simulate loudspeaker systems, which should go through most aspects of the toolbox. The first part shows how the front-end works: defining a lumped-element network, simulate the interaction between loudspeaker drivers and enclosures, setting-up exterior and interior radiation studies; the second part dives into the back-end of the toolbox and how studies can be made bypassing the main **loudspeakerSystem** class. The BEM wrapper is implemented using [bempp-cl](https://bempp.com/) which is installed alongside ElectroacPy. BEM simulations wouldn’t have been possible without the contributors of the bempp-cl project.

You'll find the GitHub repo [here](https://github.com/TchoumTchoum/electroacPy).

## Installation
### Before Starting

- The following steps have been verified on **Windows** and **Linux** for Python versions 3.9 to 3.11. For **macOS**, only version 3.9 has been tested. However, in theory, any version of Python should work as long as all dependencies are available.
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
```shell
conda create -n acoustic_sim
```

2. **Activate the environment**:
```shell
conda activate acoustic_sim
```
3. **Install Python 3.11 and pip** (you can adjust the Python version if needed):
```shell
conda install python=3.11 pip
```
4. **Install electroacPy**:

For standard installation:

```shell
pip install /path/to/electroacPy
```

For development installation:
```shell
pip install -e /path/to/electroacPy
```

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
pip install pyopencl intel-opencl-rt
```
You'll also need to install OpenCL drivers, which you'll find [here](https://www.intel.com/content/www/us/en/developer/articles/technical/intel-cpu-runtime-for-opencl-applications-with-sycl-support.html) for intel users. For more information, you can follow the **OpenCL** section from [bempp-cl installation guide](https://bempp.com/installation.html).


## Notation
In the following document, some conventions are used:

- **class**, bold statements are used for classes,
- `.function()`, inline code preceded by a dot and closed by parenthesis relates to functions,
- `.dictionary[]`, inline code preceded by a dot and closed by square brackets refers to Python dictionaries, 
- `arg`, basic inline code is for variables and function arguments.





