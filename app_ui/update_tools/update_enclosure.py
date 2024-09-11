def update_enclosure(system):
    """
    Change the frequency axis of all driver objects of a loudspeakerSystem object
    :param system:
    :param new_freq_axis:
    :return:
    """

    box_keys = ["Vb", "whichDriver", "ref2bem", "wiring", "Nd", "eta", "config"]
    box_kwargs = ["Vf", "Lp", "Sp", "rp", "Sp2", "rp2",
                  "Mmd", "Cmd", "Rmd", "Sd",
                  "Mmd2", "Cmd2", "Rmd2", "Sd2"]

    names = []
    for box in system.enclosure:
        names.append(box)

    for name in names:
        box_current = system.enclosure[name]
        tmp_attr = {}
        tmp_kwargs = {}
        for attr, value in vars(box_current).items():  # populates new attributes
            # print("{} = {}".format(attr, value))
            if attr in box_keys:
                tmp_attr[attr] = value
            elif attr in box_kwargs:
                tmp_kwargs[attr] = value

        if tmp_attr["config"] == "sealed":
            system.enclosure.pop(name)
            system.lem_enclosure_2(name, tmp_attr["Vb"], tmp_attr["eta"], tmp_attr["whichDriver"],
                                   tmp_attr["Nd"], tmp_attr["wiring"], tmp_attr["ref2bem"])

        elif tmp_attr["config"] == "vented":
            system.enclosure.pop(name)
            system.lem_enclosure_2(name, tmp_attr["Vb"], tmp_attr["eta"], tmp_attr["whichDriver"],
                                   tmp_attr["Nd"], tmp_attr["wiring"], tmp_attr["ref2bem"],
                                   Lp=tmp_kwargs["Lp"], rp=tmp_kwargs["rp"])

        elif tmp_attr["config"] == "passiveRadiator":
            system.enclosure.pop(name)
            system.lem_enclosure_2(name, tmp_attr["Vb"], tmp_attr["eta"], tmp_attr["whichDriver"],
                                   tmp_attr["Nd"], tmp_attr["wiring"], tmp_attr["ref2bem"],
                                   Mmd=tmp_kwargs["Mmd"], Rmd=tmp_kwargs["Rmd"], Cmd=tmp_kwargs["Cmd"],
                                   Sd=tmp_kwargs["Sd"])

        elif tmp_attr["config"] == "bandpass":
            system.enclosure.pop(name)
            system.lem_enclosure_2(name, tmp_attr["Vb"], tmp_attr["eta"], tmp_attr["whichDriver"],
                                   tmp_attr["Nd"], tmp_attr["wiring"], tmp_attr["ref2bem"],
                                   Lp=tmp_kwargs["Lp"], rp=tmp_kwargs["rp"], Vf=tmp_kwargs["Vf"])


        elif tmp_attr["config"] == "bandpass_2":
            system.enclosure.pop(name)
            system.lem_enclosure_2(name, tmp_attr["Vb"], tmp_attr["eta"], tmp_attr["whichDriver"],
                                   tmp_attr["Nd"], tmp_attr["wiring"], tmp_attr["ref2bem"],
                                   Lp=tmp_kwargs["Lp"], rp=tmp_kwargs["rp"], Vf=tmp_kwargs["Vf"],
                                   Lp2=tmp_kwargs["Lp2"], rp2=tmp_kwargs["rp2"])


        elif tmp_attr["config"] == "bandpass_pr":
            system.enclosure.pop(name)
            system.lem_enclosure_2(name, tmp_attr["Vb"], tmp_attr["eta"], tmp_attr["whichDriver"],
                                   tmp_attr["Nd"], tmp_attr["wiring"], tmp_attr["ref2bem"],
                                   Mmd=tmp_kwargs["Mmd"], Rmd=tmp_kwargs["Rmd"], Cmd=tmp_kwargs["Cmd"],
                                   Sd=tmp_kwargs["Sd"], Vf=tmp_kwargs["Vf"])

        elif tmp_attr["config"] == "bandpass_pr_2":
            system.enclosure.pop(name)
            system.lem_enclosure_2(name, tmp_attr["Vb"], tmp_attr["eta"], tmp_attr["whichDriver"],
                                   tmp_attr["Nd"], tmp_attr["wiring"], tmp_attr["ref2bem"],
                                   Mmd=tmp_kwargs["Mmd"], Rmd=tmp_kwargs["Rmd"], Cmd=tmp_kwargs["Cmd"],
                                   Sd=tmp_kwargs["Sd"], Vf=tmp_kwargs["Vf"],
                                   Mmd2=tmp_kwargs["Mmd2"], Rmd2=tmp_kwargs["Rmd2"], Cmd2=tmp_kwargs["Cmd2"],
                                   Sd2=tmp_kwargs["Sd2"])

