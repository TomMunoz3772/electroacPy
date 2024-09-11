def update_study(system, view_tab, info_tab, evaluation_setup):
    if bool(system.observation) is True:
        result_info_tab = info_tab.result_tab_info  # load info tab and delete / close related window
        if result_info_tab.polar_parameter_window is not None:
            result_info_tab.polar_parameter_window.close()
            result_info_tab.polar_parameter_window = None       # reset window to rebuild-it with correct frequency axis
        if result_info_tab.field_parameter_window is not None:
            result_info_tab.field_parameter_window.close()
            result_info_tab.field_parameter_window = None
        result_info_tab.listWidget.clear()
        result_info_tab.isTabInitialized = False

        system.observation.pop("GUI")
        system.acoustic_study.pop("GUI")

        ntab = view_tab.tabs.count()
        for i in range(ntab - 1, 1, -1):
            view_tab.tabs.removeTab(i)

        # reset results
        evaluation_setup["results"] = {}

        view_tab.obsView_tabs = {"FIELDS_VIEWER": {"name": [],
                                               "plotter": None}}


