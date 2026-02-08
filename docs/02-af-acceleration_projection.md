# Projection of Acceleration Data

For those having access to a laser vibrometer or other way of estimating the acceleration of surfaces, it is possible to project such data with BEM. While the method is still experimental, some fairly good results where obtained using hexagonal measurement grids. 

This chapter goes through the two methods currently available:

- import a `*.uff` data file containing both the acceleration transfer functions and points coordinates,
- manually input the measured / simulated data as two numpy arrays: transfer functions and point coordinates.


## UFF import
UFF, "universal file format", is a type of file that contains both geometry information (e.g. coordinates of measured points) and various relevant acquisition data, such as the acceleration or velocity. In electroacPy, the UFF loader will look for data that contains "Acceleration" over "Voltage" --- the transfer function of the acceleration[^noteOnUFFExport]. 

```python
system.vibrometry_data(name,        # reference for the vibrometry data
                       file_path,   # path to the acquisition data
                       rotation,    # optional, rotation of the acquisition data
                       ref2bem)     # BEM surface group reference

```

ElectroacPy achieves the projection of acquisition data by associating and interpolating the measured points to degrees of freedom on the simulated mesh. It is usually the case for acquisition data to not have the same orientation as the simulated radiating surface (for example, because of the referential of the measurement system, the driver could be positioned toward the `+z` axis, while being oriented toward `+x` in the simulation environment). Hence, a `rotation` argument allows to rotate the measured points to face the face direction as the simulated surface.

[^noteOnUFFExport]: Usually, the export options in your acquisition software should give you multiple choices in what should be saved. If not, don't hesitate to reach out. This has been tested on a Polytec laser vibrometer.



