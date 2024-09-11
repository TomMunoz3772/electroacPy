from electroacPy import load as epload
from numpy import load as npload
from os.path import join


def load_system(path, control_tab, ref_geo, ref_system, ref_evaluation):
    """
    Load folder and auto-complete all relative GUI widgets

    TODO: CAD tab, driver tab, enclosure tab, study tab

    :param ref_system:
    :param ref_geo:
    :param path:
    :param control_tab:
    :return:
    """

    system = epload(path)
    geo_npz = npload(join(path, "geometry.npz"), allow_pickle=True)
    geo = geo_npz["geo"].item()
    eval_npz = npload(join(path, "evaluation_setup.npz"), allow_pickle=True)
    evaluation = eval_npz["evaluation"].item()

    # UPDATE GEOMETRY
    fill_dictionary(ref_geo, geo)

    # UPDATE EVALUATION
    fill_dictionary(ref_evaluation, evaluation)

    # UPDATE SYSTEM
    fill_dictionary(ref_system.driver, system.driver)
    fill_dictionary(ref_system.enclosure, system.enclosure)
    fill_dictionary(ref_system.crossover, system.crossover)
    fill_dictionary(ref_system.laser_acc, system.laser_acc)
    fill_dictionary(ref_system.acoustic_study, system.acoustic_study)
    fill_dictionary(ref_system.observation, system.observation)
    fill_dictionary(ref_system.results, system.results)
    fill_dictionary(ref_system.radiator_id, system.radiator_id)

    # AUTO-COMPLETE CAD TAB
    fill_cad_tab(ref_system, ref_geo, control_tab.cad_tab)
    fill_driver_tab(ref_system, ref_geo, control_tab.driver_tab)
    fill_enclosure_tab(ref_system, ref_geo, control_tab.enclosure_tab)
    fill_study_tab(ref_evaluation, ref_geo, control_tab.study_tab)
    return None


def fill_dictionary(target_dict, source_dict):
    """
    Fills target_dict with the contents of source_dict.
    The reference of target_dict is not overwritten.
    """
    target_dict.clear()
    target_dict.update(source_dict)

def fill_cad_tab(system, geo, cad_tab):
    cad_tab.lineEdit_location.setText(geo["cad_view_path"])
    minSize = geo["minSize"]
    maxSize = geo["maxSize"]
    ppw = geo['ppw']
    factor = int(maxSize / minSize)
    freqMax = int(system.c / maxSize / ppw)
    cad_tab.lineEdit_maxFreq.setText(str(freqMax))
    cad_tab.lineEdit_ppw.setText(str(ppw))
    cad_tab.lineEdit_minFactor.setText(str(factor))

    for surface in geo["surfaceList"]:
        cad_tab.listWidget_radSurf.addItem(surface)
        cad_tab.geo["surfaceList"][surface] = geo["surfaceList"][surface]

    cad_tab.load_step_file()
    cad_tab.geo = geo

    if geo["meshPath"] is not None:
        cad_tab.mesh_geo()
    return None


def fill_driver_tab(system, geo, driver_tab):
    for driver in system.driver:
        driver_tab.listWidget.addItem(driver)



def fill_enclosure_tab(system, geo, box_tab):
    for box in system.enclosure:
        box_tab.listWidget.addItem(box)

def fill_study_tab(evaluation, geo, study_tab):
    """
    Allows to add informations into the ListWidget of study_tab
    :param evaluation:
    :param geo:
    :param study_tab:
    :return:
    """
    for eval in evaluation:
        if eval == "radiators" or eval == "boundaries" or eval == "offsets" or eval == "results":
            pass
        else:
            study_tab.list_eval.addItem(evaluation[eval]["name"])



