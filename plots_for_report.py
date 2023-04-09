import pluto
import tools
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib as mpl
import matplotlib.colors as colors
import argparse as argp
import os
import UI_helper
import Navigation_helper
from celluloid import Camera

from Global_variables import *

import TwoD_sigma
import OneD_profile
import Nbody_Characteristics
import Accretion
import migration
import Torque
import Disc_Characteristics
import Data_parser_helper

from unittest.mock import patch, create_autospec

def for_all_methods(decorator):
    def decorate(cls):
        for attr in cls.__dict__: # there's propably a better way to do this
            if callable(getattr(cls, attr)):
                setattr(cls, attr, decorator(getattr(cls, attr)))
        return cls
    return decorate


def define_new_define_save_plot_fn(prefix: str):
    def new_save_plot_fn(plots_dir: str, file_name: str, extension: str = 'png', overwrite: bool = overwrite_plots):
        save_path = '{}{}_{}.{}'.format(plots_dir, prefix, file_name, extension)
        if not overwrite_plots:
            repeated_plots = 0
            while(os.path.isfile(save_path)):
                save_path = '{}({}).{}'.format(file_name, repeated_plots, extension)
                repeated_plots += 1
        print('Saving plot in {0}'.format(save_path))
        return save_path
    return new_save_plot_fn
    

@for_all_methods(patch('plotter_helper.define_save_plot'))
@for_all_methods(patch('UI_helper.name_the_plot'))
@for_all_methods(patch('UI_helper.define_legend_name'))
@for_all_methods(patch('UI_helper.selectPlottingRange'))
@for_all_methods(patch('UI_helper.select_averaging_length'))
@for_all_methods(patch('UI_helper.selectFunctionsToRun'))
@for_all_methods(patch('UI_helper.selectObjectsToPlot'))
@for_all_methods(patch('UI_helper.selectObjectToPlot'))
@for_all_methods(patch('UI_helper.selectManyDataFilesToPlot'))
@for_all_methods(patch('UI_helper.selectDataFileToPlot'))
@for_all_methods(patch('UI_helper.select_many_data_ids_to_overlay'))
@for_all_methods(patch('UI_helper.selectDataToPlot'))
class report_plots():
    def result_3(
        mock_selectDataToPlot,
        mock_select_many_data_ids_to_overlay,
        mock_selectDataFileToPlot,
        mock_selectManyDataFilesToPlot,
        mock_selectObjectToPlot,
        mock_selectObjectsToPlot,
        mock_selectFunctionsToRun,
        mock_select_averaging_length,
        mock_selectPlottingRange,
        mock_define_legend_name,
        mock_name_the_plot,
        mock_define_save_plot
        ):
        data_ids= [
            data_id('GROUP_can_47d_gap_clear/can_47d_gap_clear_setup1', 0, 'setup 1'),
        ]
        # mock_selectDataToPlot.return_value = data_name
        mock_select_many_data_ids_to_overlay.return_value = data_ids
        # directories = Navigation_helper.Directories(data_name)
        mock_selectDataFileToPlot.return_value = 0
        mock_selectManyDataFilesToPlot.return_value = [0, 100]
        mock_selectPlottingRange.return_value = [0, 140]
        mock_selectObjectToPlot.return_value = Kep47b_astrophysical_object
        # mock_selectFunctionsToRun.return_value = 'all'
        mock_define_legend_name.return_value = [data_ids[i].legend_name for i in range(len(data_ids))]
        mock_name_the_plot.return_value = 'result_3'

        mock_define_save_plot.side_effect = define_new_define_save_plot_fn('result_3')

        mock_select_averaging_length.return_value = 250
        mock_selectDataToPlot.return_value = data_ids[0].name
        mock_selectObjectsToPlot.return_value = [2, 3], [Jupiter_astrophysical_object, Kep47b_astrophysical_object]
        Nbody_Characteristics.plot_resonance_dat()

        Nbody_Characteristics.plot_one_using_dat_planet_and_cavity()

    def result_5(
        mock_selectDataToPlot,
        mock_select_many_data_ids_to_overlay,
        mock_selectDataFileToPlot,
        mock_selectManyDataFilesToPlot,
        mock_selectObjectToPlot,
        mock_selectObjectsToPlot,
        mock_selectFunctionsToRun,
        mock_select_averaging_length,
        mock_selectPlottingRange,
        mock_define_legend_name,
        mock_name_the_plot,
        mock_define_save_plot,
        ):
        data_ids= [
            data_id('GROUP_sustain_circ/sustain_circ_setup1_final_ab', 0, 'setup 1'),
            data_id('GROUP_sustain_circ/sustain_circ_setup2_final_ab', 0, 'setup 2'),
            # data_id('GROUP_sustain_circ/sustain_circ_setup3_final_ab', 0, 'setup 3'),
            data_id('GROUP_sustain_circ/sustain_circ_setup4_final_ab', 0, 'setup 4'),
        ]
        # mock_selectDataToPlot.return_value = data_name
        mock_select_many_data_ids_to_overlay.return_value = data_ids
        # directories = Navigation_helper.Directories(data_name)
        mock_selectDataFileToPlot.return_value = 0
        mock_selectManyDataFilesToPlot.return_value = [0, 100]
        mock_selectPlottingRange.return_value = [0, 100]
        mock_selectObjectToPlot.return_value = Kep47b_astrophysical_object
        # mock_selectFunctionsToRun.return_value = 'all'
        mock_define_legend_name.return_value = [data_ids[i].legend_name for i in range(len(data_ids))]
        mock_name_the_plot.return_value = 'result_5_many_setups'

        mock_select_averaging_length.return_value = 250
        Nbody_Characteristics.plot_many_data_id_using_dat()

        data_ids.append(data_id('GROUP_planets_start_at_final_inner_orbit/no_planet_4', 0, 'no planet'))
        mock_select_averaging_length.return_value = 50
        Disc_Characteristics.plot_one_multiple_data_sets_overlayed_disc_e_avg()

        mock_select_averaging_length.return_value = 10
        Disc_Characteristics.plot_one_multiple_data_sets_gap_parameters_out()
    
    def result_6(
        mock_selectDataToPlot,
        mock_select_many_data_ids_to_overlay,
        mock_selectDataFileToPlot,
        mock_selectManyDataFilesToPlot,
        mock_selectObjectToPlot,
        mock_selectObjectsToPlot,
        mock_selectFunctionsToRun,
        mock_select_averaging_length,
        mock_selectPlottingRange,
        mock_define_legend_name,
        mock_name_the_plot,
        mock_define_save_plot,
        ):
        data_ids= [
            data_id('GROUP_sustain_circ_accrete/sustain_circ_setup3_final_ab_accrete_00069',   0, 'setup 3, high acc'),
            data_id('GROUP_sustain_circ_accrete/sustain_circ_setup3_final_ab_accrete_000069',  0, 'setup 3, mid acc'),
            data_id('GROUP_sustain_circ_accrete/sustain_circ_setup3_final_ab_accrete_0000069', 0, 'setup 3, low acc'),
            data_id('GROUP_sustain_circ_accrete/sustain_circ_setup4_final_ab_accrete_00069',   0, 'setup 4, high acc'),
            data_id('GROUP_sustain_circ_accrete/sustain_circ_setup4_final_ab_accrete_000069',  0, 'setup 4, mid acc'),
            data_id('GROUP_sustain_circ_accrete/sustain_circ_setup4_final_ab_accrete_0000069', 0, 'setup 4, low acc'),
        ]
        # mock_selectDataToPlot.return_value = data_name
        mock_select_many_data_ids_to_overlay.return_value = data_ids
        # directories = Navigation_helper.Directories(data_name)
        mock_selectDataFileToPlot.return_value = 0
        # mock_selectManyDataFilesToPlot.return_value = [0, 240]
        mock_selectPlottingRange.return_value = [0, 240]
        mock_selectObjectToPlot.return_value = Kep47b_astrophysical_object
        # mock_selectFunctionsToRun.return_value = 'all'
        mock_define_legend_name.return_value = [data_ids[i].legend_name for i in range(len(data_ids))]
        mock_name_the_plot.return_value = 'result_6_many_setups'

        mock_select_averaging_length.return_value = 200
        Nbody_Characteristics.plot_many_data_id_using_dat()


        for data in data_ids:
            mock_selectDataToPlot.return_value = data.name
            Disc_Characteristics.calculate_total_disc_accretion()

        data_ids.append(data_id('GROUP_planets_start_at_final_inner_orbit/no_planet_4', 0, 'no planet'))
        
        mock_select_averaging_length.return_value = 50
        Disc_Characteristics.plot_one_multiple_data_sets_overlayed_disc_e_avg()   

        mock_select_averaging_length.return_value = 10
        Disc_Characteristics.plot_one_multiple_data_sets_gap_parameters_out()


    def result_7(
        mock_selectDataToPlot,
        mock_select_many_data_ids_to_overlay,
        mock_selectDataFileToPlot,
        mock_selectManyDataFilesToPlot,
        mock_selectObjectToPlot,
        mock_selectObjectsToPlot,
        mock_selectFunctionsToRun,
        mock_select_averaging_length,
        mock_selectPlottingRange,
        mock_define_legend_name,
        mock_name_the_plot,
        mock_define_save_plot
        ):
        data_ids= [
            data_id('GROUP_low_alpha_low_h/control',   0, 'control'),
            # data_id('GROUP_low_alpha_low_h/less_low_h',   0, 'h=0.02'),
            data_id('GROUP_low_alpha_low_h/low_alpha',   0, r'$\alpha$=1e-5'),
            data_id('GROUP_low_alpha_low_h/low_h',   0, 'h=0.01'),
            data_id('GROUP_low_alpha_low_h/low_alpha_low_h',   0, r'$\alpha$=1e-5, h=0.01'),
        ]
        # mock_selectDataToPlot.return_value = data_name
        mock_select_many_data_ids_to_overlay.return_value = data_ids
        # directories = Navigation_helper.Directories(data_name)
        mock_selectDataFileToPlot.return_value = 0
        # mock_selectManyDataFilesToPlot.return_value = [0, 240]
        mock_selectPlottingRange.return_value = [0, 150]
        mock_selectObjectToPlot.return_value = Kep47b_astrophysical_object
        # mock_selectFunctionsToRun.return_value = 'all'
        mock_define_legend_name.return_value = [data_ids[i].legend_name for i in range(len(data_ids))]
        mock_name_the_plot.return_value = 'result_7'

        mock_select_averaging_length.return_value = 50
        Nbody_Characteristics.plot_many_data_id_using_dat()
        
        mock_select_averaging_length.return_value = 50
        Disc_Characteristics.plot_one_multiple_data_sets_overlayed_disc_e_avg()   

        mock_select_averaging_length.return_value = 10
        Disc_Characteristics.plot_one_multiple_data_sets_gap_parameters_out()



if __name__ == '__main__':
    # report_plots.result_5()
    report_plots.result_3()