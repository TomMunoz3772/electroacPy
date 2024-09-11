def update_driver(system):
    """
    Change the frequency axis of all driver objects of a loudspeakerSystem object
    :param system:
    :param new_freq_axis:
    :return:
    """

    names = []
    for driver in system.driver:
        names.append(driver)

    for name in names:
        drv_current = system.driver[name]
        U = drv_current.U
        Le = drv_current.Le
        Re = drv_current.Re
        Cms = drv_current.Cms
        Mms = drv_current.Mms
        Rms = drv_current.Rms
        Bl = drv_current.Bl
        Sd = drv_current.Sd
        ref2bem = drv_current.ref2bem
        system.driver.pop(name)
        system.lem_driver(name, U, Le, Re, Cms, Mms, Rms, Bl, Sd, ref2bem=ref2bem)
