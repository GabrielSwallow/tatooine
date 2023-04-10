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
import plotter_helper

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
                save_path = '{}{}_{}({}).{}'.format(plots_dir, prefix, file_name, repeated_plots, extension)
                repeated_plots += 1
        print('Saving plot in {0}'.format(save_path))
        return save_path
    return new_save_plot_fn

def generate_all_plots():
        report_plots.result_1()
        report_plots.result_2()
        report_plots.result_3_and_4()
        report_plots.result_5()
        report_plots.result_6()
        report_plots.result_7()
        report_plots.result_8()

data_to_plot_list = [
    Nbody_Characteristics.possible_data_to_plot.eccentricity,
    # Nbody_Characteristics.possible_data_to_plot.semi_major_axis,
]   
@for_all_methods(patch('UI_helper.select_object_config_to_plot'))
@for_all_methods(patch('UI_helper.select_a_or_e_to_plot'))
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

    def result_1(
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
        mock_select_a_or_e_to_plot,
        mock_select_object_config_to_plot,
        ):
        result_name = 'result_1'
        data_ids= [
            data_id('GROUP_type_2_migration/type_2_migration_setup1', 0, 'setup 1'),
            data_id('GROUP_type_2_migration/type_2_migration_setup2', 0, 'setup 2'),
            data_id('GROUP_type_2_migration/type_2_migration_setup3', 0, 'setup 3'),
            data_id('GROUP_type_2_migration/type_2_migration_setup4', 0, 'setup 4'),
        ]
        # mock_selectDataToPlot.return_value = data_name
        mock_select_many_data_ids_to_overlay.return_value = data_ids
        # directories = Navigation_helper.Directories(data_name)
        mock_selectDataFileToPlot.return_value = 0
        mock_selectManyDataFilesToPlot.return_value = [0, 50, 100]
        mock_selectPlottingRange.return_value = [0, 100]
        mock_selectObjectToPlot.return_value = Kep47b_astrophysical_object
        # mock_selectFunctionsToRun.return_value = 'all'
        mock_define_legend_name.return_value = [data_ids[i].legend_name for i in range(len(data_ids))]
        mock_name_the_plot.return_value = ''
        mock_define_save_plot.side_effect = define_new_define_save_plot_fn(result_name)
        mock_select_object_config_to_plot.return_value = [Kep47b_astrophysical_object, cavity_astrophysical_object]

        mock_select_averaging_length.return_value = 250
        mock_select_a_or_e_to_plot.return_value = [Nbody_Characteristics.possible_data_to_plot.semi_major_axis]
        Nbody_Characteristics.plot_many_data_id_using_dat()
        mock_select_a_or_e_to_plot.return_value = [Nbody_Characteristics.possible_data_to_plot.eccentricity]
        Nbody_Characteristics.plot_many_data_id_using_dat()
 
        # for data_to_plot in data_to_plot_list:
        #     for data in data_ids:
        #         mock_selectDataToPlot.return_value = data.name
        #         mock_define_save_plot.side_effect = define_new_define_save_plot_fn('{}_{}'.format(result_name, data.legend_name))
        #         mock_select_a_or_e_to_plot.return_value = [data_to_plot]
        #         Nbody_Characteristics.plot_one_using_dat_planet_and_cavity()

        mock_select_a_or_e_to_plot.return_value = [Nbody_Characteristics.possible_data_to_plot.semi_major_axis]
        mock_selectDataToPlot.return_value = data_ids[0].name
        mock_define_save_plot.side_effect = define_new_define_save_plot_fn('{}_{}'.format(result_name, data_ids[0].legend_name))
        Nbody_Characteristics.plot_one_using_dat_planet_and_cavity()
    
    def result_2(
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
        mock_select_a_or_e_to_plot,
        mock_select_object_config_to_plot,
        ):
        result_name = 'result_2_TEST'
        data_ids= [
            data_id('GROUP_can_47d_gap_clear/can_47d_gap_clear_setup1', 0, ''),
        ]
        mock_select_many_data_ids_to_overlay.return_value = data_ids
        directories = Navigation_helper.Directories(data_ids[0].name, save_plots_local_to_data=True)
        mock_selectDataFileToPlot.return_value = 0
        mock_selectManyDataFilesToPlot.return_value = [0, 50, 100]
        mock_selectPlottingRange.return_value = [0, 1145]
        mock_selectObjectToPlot.return_value = Kep47b_astrophysical_object
        # mock_selectFunctionsToRun.return_value = 'all'
        mock_define_legend_name.return_value = [data_ids[i].legend_name for i in range(len(data_ids))]
        mock_name_the_plot.return_value = ''
        mock_select_object_config_to_plot.return_value = [Kep47b_astrophysical_object, cavity_astrophysical_object]
        mock_selectDataToPlot.return_value = data_ids[0].name
        mock_define_save_plot.side_effect = define_new_define_save_plot_fn('{}'.format(result_name))

        mock_select_averaging_length.return_value = 250
        # mock_select_a_or_e_to_plot.return_value = [Nbody_Characteristics.possible_data_to_plot.semi_major_axis]
        # Nbody_Characteristics.plot_one_using_dat_planet_and_cavity()
        # mock_select_a_or_e_to_plot.return_value = [Nbody_Characteristics.possible_data_to_plot.eccentricity]
        # Nbody_Characteristics.plot_one_using_dat_planet_and_cavity()

        mock_selectDataFileToPlot.return_value = 0
        mock_selectManyDataFilesToPlot.return_value = [0, 500, Navigation_helper.findMaxFileNumber(directories.out_dir)]
        TwoD_sigma.plot_many()

    def result_type_2_migration(
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
        mock_select_a_or_e_to_plot,
        mock_select_object_config_to_plot,
        ):
        result_name = 'result_type_2_migration'
        data_ids= [
            data_id('GROUP_type_2_migration/type_2_migration_setup1', 0, 'setup 1'),
            data_id('GROUP_type_2_migration/type_2_migration_setup2', 0, 'setup 2'),
            data_id('GROUP_type_2_migration/type_2_migration_setup3', 0, 'setup 3'),
            data_id('GROUP_type_2_migration/type_2_migration_setup4', 0, 'setup 4'),
        ]
        # mock_selectDataToPlot.return_value = data_name
        mock_select_many_data_ids_to_overlay.return_value = data_ids
        # directories = Navigation_helper.Directories(data_name)
        mock_selectDataFileToPlot.return_value = 0
        mock_selectManyDataFilesToPlot.return_value = [0, 50, 100]
        mock_selectPlottingRange.return_value = [0, 100]
        mock_selectObjectToPlot.return_value = Kep47b_astrophysical_object
        # mock_selectFunctionsToRun.return_value = 'all'
        mock_define_legend_name.return_value = [data_ids[i].legend_name for i in range(len(data_ids))]
        mock_name_the_plot.return_value = ''
        mock_define_save_plot.side_effect = define_new_define_save_plot_fn(result_name)
        mock_select_a_or_e_to_plot.return_value = [Nbody_Characteristics.possible_data_to_plot.semi_major_axis]
        mock_select_object_config_to_plot.return_value = [Kep47b_astrophysical_object, cavity_astrophysical_object]

        mock_select_averaging_length.return_value = 250
        mock_select_a_or_e_to_plot.return_value = [Nbody_Characteristics.possible_data_to_plot.semi_major_axis]
        Nbody_Characteristics.plot_many_data_id_using_dat()
        mock_select_a_or_e_to_plot.return_value = [Nbody_Characteristics.possible_data_to_plot.eccentricity]
        Nbody_Characteristics.plot_many_data_id_using_dat()
 
        for data_to_plot in data_to_plot_list:
            for data in data_ids:
                mock_selectDataToPlot.return_value = data.name
                mock_define_save_plot.side_effect = define_new_define_save_plot_fn('{}_{}'.format(result_name, data.legend_name))
                mock_select_a_or_e_to_plot.return_value = [data_to_plot]
                Nbody_Characteristics.plot_one_using_dat_planet_and_cavity()


        mock_select_a_or_e_to_plot.return_value = [Nbody_Characteristics.possible_data_to_plot.semi_major_axis]
        mock_selectDataToPlot.return_value = data_ids[0].name
        mock_define_save_plot.side_effect = define_new_define_save_plot_fn('{}_{}'.format(result_name, data_ids[0].legend_name))
        Nbody_Characteristics.plot_one_using_dat_planet_and_cavity()

    def result_3_and_4(
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
        mock_select_a_or_e_to_plot,
        mock_select_object_config_to_plot,
        ):
        result_name = 'result_3_and_4'
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
        mock_name_the_plot.return_value = result_name
        mock_define_save_plot.side_effect = define_new_define_save_plot_fn(result_name)
        mock_select_a_or_e_to_plot.return_value = [Nbody_Characteristics.possible_data_to_plot.semi_major_axis]

        Jupiter_astrophysical_object.id = 2
        Kep47b_astrophysical_object.id = 3
        mock_select_object_config_to_plot.return_value = [Jupiter_astrophysical_object, Kep47b_astrophysical_object] 

        mock_select_averaging_length.return_value = 250
        mock_selectDataToPlot.return_value = data_ids[0].name
        mock_selectObjectsToPlot.return_value = [2, 3], [Jupiter_astrophysical_object, Kep47b_astrophysical_object]
        Nbody_Characteristics.plot_resonance_dat()

        mock_select_a_or_e_to_plot.return_value = [Nbody_Characteristics.possible_data_to_plot.semi_major_axis]
        Nbody_Characteristics.plot_one_using_dat_planet_and_cavity()

        Jupiter_astrophysical_object.id = 7
        Kep47b_astrophysical_object.id = 2

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
        mock_select_a_or_e_to_plot,
        mock_select_object_config_to_plot,
        ):
        result_name = 'result_5'
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
        mock_name_the_plot.return_value = ''
        mock_define_save_plot.side_effect = define_new_define_save_plot_fn(result_name)
        mock_select_a_or_e_to_plot.return_value = [Nbody_Characteristics.possible_data_to_plot.semi_major_axis]
        mock_select_object_config_to_plot.return_value = [Kep47b_astrophysical_object] 

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
        mock_select_a_or_e_to_plot,
        mock_select_object_config_to_plot,
        ):
        result_name = 'result_6'

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
        mock_name_the_plot.return_value = ''
        mock_define_save_plot.side_effect = define_new_define_save_plot_fn(result_name)
        mock_select_a_or_e_to_plot.return_value = [Nbody_Characteristics.possible_data_to_plot.semi_major_axis]
        mock_select_object_config_to_plot.return_value = [Kep47b_astrophysical_object] 

        mock_select_averaging_length.return_value = 200
        Nbody_Characteristics.plot_many_data_id_using_dat()

        with open('{}{}_accretion_data.txt'.format(global_plots_dir, result_name), 'w') as f:
            for data in data_ids:
                mock_selectDataToPlot.return_value = data.name
                acc_tot, print_statement = Disc_Characteristics.calculate_total_disc_accretion()
                f.write(print_statement + '\n')

        # data_ids.append(data_id('GROUP_planets_start_at_final_inner_orbit/no_planet_4', 0, 'no planet'))
        
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
        mock_define_save_plot,
        mock_select_a_or_e_to_plot,
        mock_select_object_config_to_plot,
        ):
        result_name = 'result_7'

        data_ids= [
            data_id('GROUP_can_47d_keep_47b_stable_BETTER/setup1',   0, 'setup1'),
            data_id('GROUP_can_47d_keep_47b_stable_BETTER/setup2',   0, 'setup2'),
            data_id('GROUP_can_47d_keep_47b_stable_BETTER/setup3',   0, 'setup3'),
            data_id('GROUP_can_47d_keep_47b_stable_BETTER/setup4',   0, 'setup4'),
        ]
        # mock_selectDataToPlot.return_value = data_name
        mock_select_many_data_ids_to_overlay.return_value = data_ids
        # directories = Navigation_helper.Directories(data_name)
        mock_selectDataFileToPlot.return_value = 0
        # mock_selectManyDataFilesToPlot.return_value = [0, 240]
        mock_selectPlottingRange.return_value = [0, 100]
        # mock_selectFunctionsToRun.return_value = 'all'
        mock_define_legend_name.return_value = [data_ids[i].legend_name for i in range(len(data_ids))]
        mock_name_the_plot.return_value = ''
        mock_define_save_plot.side_effect = define_new_define_save_plot_fn(result_name)
        mock_select_a_or_e_to_plot.return_value = [Nbody_Characteristics.possible_data_to_plot.semi_major_axis]
        mock_select_object_config_to_plot.return_value = [Kep47b_astrophysical_object, Kep47d_astrophysical_object] 

        mock_selectObjectToPlot.return_value = Kep47b_astrophysical_object
        mock_select_averaging_length.return_value = 50
        Nbody_Characteristics.plot_many_data_id_using_dat()

        mock_selectObjectToPlot.return_value = Kep47d_astrophysical_object
        mock_select_averaging_length.return_value = 50
        Nbody_Characteristics.plot_many_data_id_using_dat()
        
        # mock_select_averaging_length.return_value = 50
        # Disc_Characteristics.plot_one_multiple_data_sets_overlayed_disc_e_avg()   

        # mock_select_averaging_length.return_value = 10
        # Disc_Characteristics.plot_one_multiple_data_sets_gap_parameters_out()

    def result_8(
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
        mock_select_a_or_e_to_plot,
        mock_select_object_config_to_plot,
        ):
        result_name = 'result_8'

        data_ids= [
            data_id('GROUP_low_alpha_low_h/control',   0, 'control'),
            data_id('GROUP_low_alpha_low_h/less_low_h',   0, 'h=0.02'),
            data_id('GROUP_low_alpha_low_h/low_alpha',   0, r'$\alpha$=1e-5'),
            data_id('GROUP_low_alpha_low_h/low_h',   0, 'h=0.01'),
            data_id('GROUP_low_alpha_low_h/low_alpha_low_h',   0, r'$\alpha$=1e-5, h=0.01'),
        ]
        # mock_selectDataToPlot.return_value = data_name
        mock_select_many_data_ids_to_overlay.return_value = data_ids
        # directories = Navigation_helper.Directories(data_name)
        mock_selectDataFileToPlot.return_value = 0
        # mock_selectManyDataFilesToPlot.return_value = [0, 240]
        mock_selectPlottingRange.return_value = [0, 111]
        mock_selectObjectToPlot.return_value = Kep47b_astrophysical_object
        # mock_selectFunctionsToRun.return_value = 'all'
        mock_define_legend_name.return_value = [data_ids[i].legend_name for i in range(len(data_ids))]
        mock_name_the_plot.return_value = ''
        mock_define_save_plot.side_effect = define_new_define_save_plot_fn(result_name)
        mock_select_a_or_e_to_plot.return_value = [Nbody_Characteristics.possible_data_to_plot.semi_major_axis]
        mock_select_object_config_to_plot.return_value = [Kep47b_astrophysical_object] 

        mock_select_averaging_length.return_value = 50
        Nbody_Characteristics.plot_many_data_id_using_dat()
        
        mock_select_averaging_length.return_value = 50
        Disc_Characteristics.plot_one_multiple_data_sets_overlayed_disc_e_avg()   

        mock_select_averaging_length.return_value = 10
        Disc_Characteristics.plot_one_multiple_data_sets_gap_parameters_out()
        

if __name__ == '__main__':
    # report_plots.result_3_and_4()
    report_plots.result_2()
    # generate_all_plots()