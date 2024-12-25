# LEM modules

## speakerSim

In electroacPy, all speaker-related classes are grouped within the **speakerSim** submodule:

- **electroAcousticDriver**: basic class to store driver parameters and plot related impedance and mechanical behavior --- where functions like `.plotZe()` and `.plotXVA()` are defined. These objects are stored in the `.driver[]` dictionary of loudspeaker systems,
- **speakerBox**: class to create enclosures based on given dimensions. Can be associated to drivers to estimate related velocities. **speakerBox** are saved in the `.enclosure[]` dictionary.
- **xovers**: filtering functions applied to radiated pressure. Stored in the `.network[]` dictionary.

### **electroAcousticDriver**

When creating a driver object, the circuit shown in {numref}`ead-model` is used to compute its free-air electrical impedance from electrical and mechanical behavior[^driverLimit]. The object itself receives the `EAC` identifier for the **loudspeakerSystem** class.

[^driverLimit]: For now the only lumped element model supported in electroacPy is [Le, Re], if necessary, other models will be added in the future.


```{figure} ./LEM_backend_images/00_electroacoustic_driver.svg
    :name: ead-model

Lumped Element Model used in **electroAcousticDriver**.
```

From this initial circuit we get

$$
Z_{in} = Z_e + \frac{(Bl)^2}{Z_{ms}},
$$

with $Z_e$ the coil impedance expressed as

$$
Z_e = R_e + j\omega L_e
$$

and $Z_{ms}$ the mechanical impedance of moving parts: 

$$
Z_{ms} = R_{ms} + j\omega M_{ms} + \frac{1}{j\omega C_{ms}}.
$$

A good thing to note is that `.plotZe()` will actually plot the input impedance $Z_{in}$ and not the coil impedance $Z_e$.


### **speakerBox**
As explained in the lumped-element modeler for enclosures, all **speakerBox** objects are primarily defined by their respective box volumes. Then, depending on user inputs, different types of systems can be created. This section illustrates the various lumped networks for each enclosure design and clarifies which parameters correspond to specific elements of the system.

#### sealed enclosure
The sealed box is the default configuration. In this case, the back of the driver is loaded by a simple acoustic volume.

```{figure} ./LEM_backend_images/0a_sealed_box.svg
    :name: sealed-box-a
    :scale: 50%  

Physical representation of the sealed box.
```

```python
system.lem_enclosure("sealed", Vb)
```


```{figure} ./LEM_backend_images/01_sealed_enclosure.svg
    :name: sealed-box

Acoustic network of a sealed enclosure. Connector `C` is linked to the front radiation impedance, while `D` is connected to the back radiation of the driver.
```

This acoustic volume is represented as a capacitance $C_{ab}$​, coupled with two resistances, $R_{ab}$​ and $R_{al}$​, which account for damping within the enclosure and acoustic leakage, respectively. These elements are defined as:

$$

C_{ab} = \frac{V_b}{\rho c}, ~~ R_{ab} = \frac{1}{2 \pi f_c C_{ab} Q_{ab}}, ~~ R_{al} = \frac{Q_l}{2\pi f_c C_{ab}},

$$

with:

- $V_b$ the enclosure volume, 
- $rho$, $c$ the air density and speed of sound,
- $f_c$ the resonance frequency of the driver inside the enclosure, 
- $Q_{ab}$ the quality factor linked to acoustic losses in the enclosure,
- $Q_{al}$ the quality factor linked to acoustic leakage out of the enclosure.

By default in the `.lem_enclosure()` method, the Q factors are set to be `Qab=120` and `Qal=30`. These values usually sits between:

$$
30<Q_{ab}<100, ~~ 10<Q_{al}<30,
$$

where lower values on the left indicate higher losses, while higher values on the right correspond to lower losses. Of course, these values depend heavily on the amount of damping material in the enclosure and the quality of construction. You will surely need to determine the "correct" values through trial and error, depending on your build.

#### ported enclosure
Ported enclosures[^ported-box] are defined by adding --- to the `.lem_enclosure()` method --- a port's length $L_p$ and a corresponding radius *or* cross-section area: $r_p$ or $S_p$, respectively. Since ports are more than just straight tubes, it is possible to specify whether none, one, or both ends are flanged. This can be adjusted by setting the `flange` argument to either `"none"`, `"single"` or `"both"`. The following code shows the basics arguments for defining a ported box:

[^ported-box]: Also known as bass-reflex.

```{figure} ./LEM_backend_images/0b_ported_box.svg
    :name: ported-box-b
    :scale: 50%

Representation of a bass-reflex system with a port at the back of the enclosure.
```

```python
system.lem_enclosure("port", Vb, Lp, rp, flange="single")
```

```{figure} ./LEM_backend_images/02_ported_enclosure.svg
    :name:ported-box

Ported enclosure. The port radiation is modeled with an additional radiation impedance.
```

Relative to the acoustic network, a ported enclosure adds an acoustic mass and resistance linked to the port dimensions and another radiation impedance $Z_{radp}$ linked to the area of the opening. These are defined as:

$$
M_p = \frac{\rho}{Sp} \left( Lp + \sigma \sqrt{Sp} \right) 
$$


$$
R_p = \frac{\sqrt{2\omega \rho \mu}}{Sp} \times \left( \frac{Lp}{rp} + 1\right)
$$

where 

- $\rho$ is the air density,
- $Lp$, $rp$ and $Sp$ are the length, radius and cross-sectional area of the port, 
- $\sigma$ is the flange coefficient: 0.96 for both ends flanged, 0.84 for one flanged end, and 0.72 for both ends unflanged,
- $\omega$ is the angular frequency, 
- $\mu$ is the air viscosity (approximately $1.82\times 10^{-5}$ kg/(m.s) at 20°C).



#### ABR configuration

The ABR (or passive radiator) setup is really similar to the ported enclosure. The only difference is that it replaces the port (acoustic mass and resistance) by a membrane. Hence, the ABR is expressed as mechanical mass, resistance and compliance (`Mmd`, `Rmd`, `Cmd` and associated radiating surface `Sd`) --- similar to what is found in the mechanical circuit of a electro-dynamic loudspeaker.

```{figure} ./LEM_backend_images/0c_abr_box.svg
    :name: abr-box-c
    :scale: 50%

ABR setup, with the auxiliary radiator represented in red.
```

```python
system.lem_enclosure("abr", Vb, Mmd, Cmd, Rmd, Sd)
```

```{figure} ./LEM_backend_images/03_abr_enclosure.svg
    :name:abr-box

ABR enclosure. Here, the port is replaced by a mechanical radiator expressed in the acoustic domain.
```

These mechanical parameters are given by:

$$

Mma = \frac{Mmd}{Sd^2}, ~~ Rma = \frac{Rmd}{Sd^2}, ~~ Cma = Cmd \times Sd^2,

$$

where $Mmd$ is the mass of the passive radiator in $kg$, $Rmd$ the mechanical resistance of suspensions in $N.s/m$ and Cmd is the suspensions compliance in $m/N$. Finally, $Sd$ is the radiating surface of the ABR, which is used to translate the mechanical elements into their acoustical counterparts.


#### 4th order bandpass (port / ABR)
The 4th order bandpass enclosure has its loudspeaker radiating into a front volume $V_f$. In that configuration, the port (or ABR) is the only radiating element of the system --- connected to the front volume. Similarly to bass-reflex or ABR enclosures, the port dimensions or ABR parameters must be given as key arguments.


```{figure} ./LEM_backend_images/0d_bandpassp_box.svg
    :name: bp-box-d
    :scale: 50%
    :align: left

Ported bandpass enclosure.
```

```{figure} ./LEM_backend_images/0e_bandpassabr_box.svg
    :name: abr-box-d
    :scale: 50%
    :align: right

ABR bandpass enclosure.
```
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>

```python
system.lem_enclosure("bp_port", Vb, Vf, Lp, rp)
system.lem_enclosure("bp_abr", Vb, Vf, Mmd, Cmd, Rmd, Sd)
```

```{figure} ./LEM_backend_images/04_bandpass_port.svg
    :name:bp-box

4th order bandpass enclosure. The port is in the volume in front of the woofer.
```

```{figure} ./LEM_backend_images/05_bandpass_abr.svg
    :name:bp-abr-box

4th order bandpass enclosure. A passive radiator replaces the port.
```


#### 6th order bandpass (port / ABR)


```{figure} ./LEM_backend_images/0h_bandpassp6_box.svg
    :name: bp6-box-h
    :scale: 50%
    :align: left

6th order bandpass --- port.
```

```{figure} ./LEM_backend_images/0g_bandpassabr6_box.svg
    :name: abr6-box-g
    :scale: 50%
    :align: right

6th order bandpass --- ABR.
```
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>

```python
system.lem_enclosure("bp_port", Vb, Vf, Lp, rp, 
                                        Lp2, rp2)
system.lem_enclosure("bp_abr", Vb, Vf, Mmd, Cmd, Rmd, Sd, 
                                       Mmd2, Cmd2, Rmd2, Sd2)
```

```{figure} ./LEM_backend_images/06_bandpass6_port.svg
    :name:bp6-port-box
    :align: center

6th order bandpass --- port.
```

```{figure} ./LEM_backend_images/07_bandpass6_abr.svg
    :name:bp6-abr-box
    :align: center

6th order bandpass --- ABR.
```
