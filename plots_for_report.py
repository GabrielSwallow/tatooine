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
    # report_plots.result_1()
    # report_plots.result_type_2_migration()
    # report_plots.result_3()
    report_plots.result_4()
    report_plots.result_5()
    report_plots.result_6a()
    report_plots.result_6b()
    report_plots.result_7()

def generate_all_LUKE_plots():
    report_plots.result_1_LUKE()
    report_plots.result_2_LUKE()
    report_plots.result_3_LUKE()

data_to_plot_list = [
    Nbody_Characteristics.possible_data_to_plot.eccentricity,
    Nbody_Characteristics.possible_data_to_plot.semi_major_axis,
]   
@for_all_methods(patch('UI_helper.show_instability_limit'))
@for_all_methods(patch('UI_helper.show_47b_final_orbit'))
@for_all_methods(patch('UI_helper.define_size_of_plot_in_abin'))
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

    def intro_disc_distribution(mock_selectDataToPlot,
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
        mock_define_size_of_plot_in_abin,
        mock_show_47b_final_orbit,
        mock_show_instability_limit,
        ):
        plot_name = 'intro_disc_distribution'
        data = data_id('GROUP_init_disc/init_disc_setup4', 0, 'setup 4')
        mock_selectDataToPlot.return_value = data.name
        mock_selectDataFileToPlot.return_value = 0
        mock_define_size_of_plot_in_abin.return_value = 12.
        mock_define_save_plot.side_effect = define_new_define_save_plot_fn('{}_{}'.format(plot_name, data.legend_name))
        TwoD_sigma.plot_one()


    def theory(
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
        mock_define_size_of_plot_in_abin,
        mock_show_47b_final_orbit,
        mock_show_instability_limit,
        ):
        plot_name = 'theory'
        data_ids= [
            data_id('GROUP_init_disc/init_disc_setup4', 0, 'init_disc_s4'),
        ]
        mock_select_many_data_ids_to_overlay.return_value = data_ids
        mock_selectDataFileToPlot.return_value = 0
        mock_selectDataToPlot.return_value = data_ids[0].name
        mock_define_save_plot.side_effect = define_new_define_save_plot_fn('{}_{}'.format(plot_name, data_ids[0].legend_name))
        mock_define_size_of_plot_in_abin.return_value = 40.
        # Disc_Characteristics.plot_difference_between_min_and_max_rho_same_radius()

        data_ids= [
            data_id('GROUP_types_of_migration/Earth', 0, 'Earth'),
            data_id('GROUP_types_of_migration/Jupiter', 0, 'Jupiter'),
            data_id('GROUP_types_of_migration/Neptune', 0, 'Neptune'),
        ]
        mock_define_size_of_plot_in_abin.return_value = 18.
        for dat in data_ids:
            mock_selectDataToPlot.return_value = dat.name
            mock_define_save_plot.side_effect = define_new_define_save_plot_fn('{}_{}'.format(plot_name, dat.legend_name))
            directories = Navigation_helper.Directories(dat.name, save_plots_local_to_data=True)
            mock_selectDataFileToPlot.return_value = Navigation_helper.findMaxFileNumber(directories.out_dir)
            TwoD_sigma.plot_one()


    def result_1_LUKE(
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
        mock_define_size_of_plot_in_abin,
        mock_show_47b_final_orbit,
        mock_show_instability_limit,
        ):
        result_name = 'LUKE_result_1'
        data_ids= [
            data_id('GROUP_disc_alp1e-3_tests/141122_1_restart', 0, ''),
        ]
        mock_select_many_data_ids_to_overlay.return_value = data_ids
        directories = Navigation_helper.Directories(data_ids[0].name, save_plots_local_to_data=True)
        mock_selectDataFileToPlot.return_value = 0
        mock_selectPlottingRange.return_value = [0, Navigation_helper.findMaxFileNumber(directories.out_dir)]

        heavy_planet_object = astrophysical_object(2, 'heavy type I', 'heavy type I', 0.)

        mock_selectObjectToPlot.return_value = heavy_planet_object
        mock_define_legend_name.return_value = [data_ids[i].legend_name for i in range(len(data_ids))]
        mock_name_the_plot.return_value = ''
        mock_select_object_config_to_plot.return_value = [heavy_planet_object, cavity_astrophysical_object]
        mock_select_averaging_length.return_value = 250
        mock_define_size_of_plot_in_abin.return_value = 12.

        mock_show_47b_final_orbit.return_value = True
        mock_show_instability_limit.return_value = False

        for data_to_plot in data_to_plot_list:
            mock_selectDataToPlot.return_value = data_ids[0].name
            mock_define_save_plot.side_effect = define_new_define_save_plot_fn('{}'.format(result_name))
            mock_select_a_or_e_to_plot.return_value = [data_to_plot]
            Nbody_Characteristics.plot_one_using_dat_planet_and_cavity()
        
        mock_selectManyDataFilesToPlot.return_value = [0, Navigation_helper.findMaxFileNumber(directories.out_dir)]
        TwoD_sigma.plot_many()
    
    def result_2_LUKE(
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
        mock_define_size_of_plot_in_abin,
        mock_show_47b_final_orbit,
        mock_show_instability_limit,
        ):
        result_name = 'LUKE_result_2'
        data_ids= [
            data_id('GROUP_disc_alp1e-3_tests/MJ1.0_aB13.0_Ve-3_restart5', 0, ''),
        ]
        mock_select_many_data_ids_to_overlay.return_value = data_ids
        directories = Navigation_helper.Directories(data_ids[0].name, save_plots_local_to_data=True)
        mock_selectDataFileToPlot.return_value = 0

        Jupiter_astrophysical_object.id = 2
        mock_selectObjectToPlot.return_value = Jupiter_astrophysical_object

        # mock_selectFunctionsToRun.return_value = 'all'
        mock_define_legend_name.return_value = [data_ids[i].legend_name for i in range(len(data_ids))]
        mock_name_the_plot.return_value = ''
        mock_select_object_config_to_plot.return_value = [Jupiter_astrophysical_object, cavity_astrophysical_object]
        mock_selectDataToPlot.return_value = data_ids[0].name
        mock_define_save_plot.side_effect = define_new_define_save_plot_fn('{}'.format(result_name))

        mock_show_47b_final_orbit.return_value = False
        mock_show_instability_limit.return_value = True
        mock_define_size_of_plot_in_abin.return_value = 12.

        ###
        # Don't want to plot max in Luke's data, as it includes ejection
        mock_selectPlottingRange.return_value = [0, 1145]
        mock_selectDataFileToPlot.return_value = 1145
        ###
        mock_select_averaging_length.return_value = 250
        mock_select_a_or_e_to_plot.return_value = [Nbody_Characteristics.possible_data_to_plot.semi_major_axis]
        Nbody_Characteristics.plot_one_using_dat_planet_and_cavity()
        mock_select_a_or_e_to_plot.return_value = [Nbody_Characteristics.possible_data_to_plot.eccentricity]
        Nbody_Characteristics.plot_one_using_dat_planet_and_cavity()


        TwoD_sigma.plot_one()

        mock_define_size_of_plot_in_abin.return_value = 12.
        Disc_Characteristics.plot_difference_between_min_and_max_rho_same_radius()

        Jupiter_astrophysical_object.id = 7 

    def result_3_LUKE(
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
        mock_define_size_of_plot_in_abin,
        mock_show_47b_final_orbit,
        mock_show_instability_limit,
        ):
        result_name = 'LUKE_result_3'
        data_ids= [
            data_id('GROUP_disc_alp1e-3_tests/GROUP_circularized_disc_tests/281122_MJ1.0_MJ0.1', 0, ''),
        ]
        mock_select_many_data_ids_to_overlay.return_value = data_ids
        mock_selectDataToPlot.return_value = data_ids[0].name
        directories = Navigation_helper.Directories(data_ids[0].name, save_plots_local_to_data=True)
        mock_selectDataFileToPlot.return_value = 0
        mock_selectManyDataFilesToPlot.return_value = [0, 100]
        mock_selectObjectToPlot.return_value = Kep47b_astrophysical_object
        mock_define_legend_name.return_value = [data_ids[i].legend_name for i in range(len(data_ids))]
        mock_name_the_plot.return_value = result_name
        mock_define_save_plot.side_effect = define_new_define_save_plot_fn(result_name)
        mock_select_a_or_e_to_plot.return_value = [Nbody_Characteristics.possible_data_to_plot.semi_major_axis]

        Jupiter_astrophysical_object.id = 2
        Kep47b_astrophysical_object.id = 3
        mock_select_object_config_to_plot.return_value = [Jupiter_astrophysical_object, Kep47b_astrophysical_object] 
        mock_selectObjectsToPlot.return_value = [2, 3], [Jupiter_astrophysical_object, Kep47b_astrophysical_object]

        mock_show_47b_final_orbit.return_value = False
        mock_show_instability_limit.return_value = True
        mock_define_size_of_plot_in_abin.return_value = 12.

        ###
        # Don't want to plot max in Luke's data, as it includes ejection
        mock_selectPlottingRange.return_value = [0, 665]
        ###
        mock_select_averaging_length.return_value = 250
        Nbody_Characteristics.plot_resonance_dat()
        mock_select_a_or_e_to_plot.return_value = [Nbody_Characteristics.possible_data_to_plot.semi_major_axis]
        Nbody_Characteristics.plot_one_using_dat_planet_and_cavity()

        Jupiter_astrophysical_object.id = 7
        Kep47b_astrophysical_object.id = 2


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
        mock_define_size_of_plot_in_abin,
        mock_show_47b_final_orbit,
        mock_show_instability_limit,
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
        mock_define_size_of_plot_in_abin,
        mock_show_47b_final_orbit,
        mock_show_instability_limit,
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
        mock_define_save_plot,
        mock_select_a_or_e_to_plot,
        mock_select_object_config_to_plot,
        mock_define_size_of_plot_in_abin,
        mock_show_47b_final_orbit,
        mock_show_instability_limit,
        ):
        result_name = 'result_3'
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

    def result_4(
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
        mock_define_size_of_plot_in_abin,
        mock_show_47b_final_orbit,
        mock_show_instability_limit,
        ):
        result_name = 'result_4_FULL'
        data_ids= [
            data_id('GROUP_sustain_circ_FULL/sustain_circ_setup1_final_ab', 0, 'setup 1'),
            data_id('GROUP_sustain_circ_FULL/sustain_circ_setup2_final_ab', 0, 'setup 2'),
            data_id('GROUP_sustain_circ_FULL/sustain_circ_setup3_final_ab', 0, 'setup 3'),
            data_id('GROUP_sustain_circ_FULL/sustain_circ_setup4_final_ab', 0, 'setup 4'),
        ]
        # mock_selectDataToPlot.return_value = data_name
        mock_select_many_data_ids_to_overlay.return_value = data_ids
        # directories = Navigation_helper.Directories(data_name)
        mock_selectDataFileToPlot.return_value = 0
        mock_selectManyDataFilesToPlot.return_value = [0, 100]
        mock_selectPlottingRange.return_value = [0, 249]
        mock_selectObjectToPlot.return_value = Kep47b_astrophysical_object
        # mock_selectFunctionsToRun.return_value = 'all'
        mock_define_legend_name.return_value = [data_ids[i].legend_name for i in range(len(data_ids))]
        mock_name_the_plot.return_value = ''
        mock_define_save_plot.side_effect = define_new_define_save_plot_fn(result_name)
        mock_select_a_or_e_to_plot.return_value = [Nbody_Characteristics.possible_data_to_plot.semi_major_axis]
        mock_select_object_config_to_plot.return_value = [Kep47b_astrophysical_object] 

        mock_show_47b_final_orbit.return_value = False
        mock_show_instability_limit.return_value = False

        mock_select_averaging_length.return_value = 250
        Nbody_Characteristics.plot_many_data_id_using_dat()

        # data_ids.append(data_id('GROUP_planets_start_at_final_inner_orbit/no_planet_4', 0, 'setup 4, no planet'))
        mock_select_averaging_length.return_value = 50
        Disc_Characteristics.plot_one_multiple_data_sets_overlayed_disc_e_avg()

        mock_select_averaging_length.return_value = 10
        Disc_Characteristics.plot_one_multiple_data_sets_gap_parameters_out()
    
    def result_4b(
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
        mock_define_size_of_plot_in_abin,
        mock_show_47b_final_orbit,
        mock_show_instability_limit,
        ):
        result_name = 'result_4b'
        data_ids= [
            data_id('GROUP_SC_large_mass_FULL/setup2', 0, r'$40 M_{\oplus}$, setup 2'), #, 'tab:blue'),
            data_id('GROUP_SC_large_mass_FULL/setup4', 0, r'$40 M_{\oplus}$, setup 4'), #, 'tab:orange'),
            data_id('GROUP_SC_large_mass_FULL/setup1', 0, r'$40 M_{\oplus}$, setup 1'), #, 'tab:green'),
            data_id('GROUP_SC_large_mass_FULL/setup3', 0, r'$40 M_{\oplus}$, setup 3'), #, 'tab:red'),
            data_id('GROUP_sustain_circ_FULL/sustain_circ_setup2_final_ab', 0, r'$2 M_{\oplus}$, setup 2'), #, 'tab:blue'),
            data_id('GROUP_sustain_circ_FULL/sustain_circ_setup4_final_ab', 0, r'$2 M_{\oplus}$, setup 4'), #, 'tab:orange'),
            data_id('GROUP_sustain_circ_FULL/sustain_circ_setup1_final_ab', 0, r'$2 M_{\oplus}$, setup 1'), #, 'tab:green'),
            data_id('GROUP_sustain_circ_FULL/sustain_circ_setup3_final_ab', 0, r'$2 M_{\oplus}$, setup 3'), #, 'tab:red'),
        ]
        # mock_selectDataToPlot.return_value = data_name
        mock_select_many_data_ids_to_overlay.return_value = data_ids
        # directories = Navigation_helper.Directories(data_name)
        mock_selectDataFileToPlot.return_value = 0
        mock_selectManyDataFilesToPlot.return_value = [0, 100]
        mock_selectPlottingRange.return_value = [0, 206]
        mock_selectObjectToPlot.return_value = Kep47b_astrophysical_object
        # mock_selectFunctionsToRun.return_value = 'all'
        mock_define_legend_name.return_value = [data_ids[i].legend_name for i in range(len(data_ids))]
        mock_name_the_plot.return_value = ''
        mock_define_save_plot.side_effect = define_new_define_save_plot_fn(result_name)
        mock_select_a_or_e_to_plot.return_value = [Nbody_Characteristics.possible_data_to_plot.semi_major_axis]
        mock_select_object_config_to_plot.return_value = [Kep47b_astrophysical_object] 

        mock_show_47b_final_orbit.return_value = False
        mock_show_instability_limit.return_value = False
        mock_define_size_of_plot_in_abin.return_value = 15.

        # mock_select_averaging_length.return_value = 250
        # Nbody_Characteristics.plot_many_data_id_using_dat()

        # # data_ids.append(data_id('GROUP_planets_start_at_final_inner_orbit/no_planet_4', 0, 'setup 4, no planet'))
        # # mock_select_averaging_length.return_value = 50
        # # Disc_Characteristics.plot_one_multiple_data_sets_overlayed_disc_e_avg()
        # # change order of plot for legend
        # data_ids= [
        #     data_id('GROUP_SC_large_mass_FULL/setup2', 0, r'$40 M_{\oplus}$, setup 2'), #, 'tab:blue'),
        #     data_id('GROUP_sustain_circ_FULL/sustain_circ_setup2_final_ab', 0, r'$2 M_{\oplus}$, setup 2'), #, 'tab:blue'),
        #     data_id('GROUP_sustain_circ_FULL/sustain_circ_setup4_final_ab', 0, r'$2 M_{\oplus}$, setup 4'), #, 'tab:orange'),
        #     data_id('GROUP_SC_large_mass_FULL/setup4', 0, r'$40 M_{\oplus}$, setup 4'), #, 'tab:orange'),
        #     data_id('GROUP_SC_large_mass_FULL/setup1', 0, r'$40 M_{\oplus}$, setup 1'), #, 'tab:green'),
        #     data_id('GROUP_sustain_circ_FULL/sustain_circ_setup1_final_ab', 0, r'$2 M_{\oplus}$, setup 1'), #, 'tab:green'),
        #     data_id('GROUP_SC_large_mass_FULL/setup3', 0, r'$40 M_{\oplus}$, setup 3'), #, 'tab:red'),
        #     data_id('GROUP_sustain_circ_FULL/sustain_circ_setup3_final_ab', 0, r'$2 M_{\oplus}$, setup 3'), #, 'tab:red'),
        # ]
        # mock_select_many_data_ids_to_overlay.return_value = data_ids

        # mock_select_a_or_e_to_plot.return_value = [Nbody_Characteristics.possible_data_to_plot.eccentricity]
        # mock_select_averaging_length.return_value = 10
        # Disc_Characteristics.plot_one_multiple_data_sets_gap_parameters_out()

        mock_define_save_plot.side_effect = define_new_define_save_plot_fn('{}_large_mass'.format(result_name))
        mock_selectDataFileToPlot.return_value = 250
        mock_selectDataToPlot.return_value = 'GROUP_SC_large_mass_FULL/setup2'
        TwoD_sigma.plot_one()

        mock_define_save_plot.side_effect = define_new_define_save_plot_fn('{}_small_mass'.format(result_name))
        mock_selectDataFileToPlot.return_value = 251
        mock_selectDataToPlot.return_value = 'GROUP_sustain_circ_FULL/sustain_circ_setup2_final_ab'
        TwoD_sigma.plot_one()
    
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
        mock_define_size_of_plot_in_abin,
        mock_show_47b_final_orbit,
        mock_show_instability_limit,
        ):
        result_name = 'result_5'

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

        mock_show_47b_final_orbit.return_value = False
        mock_show_instability_limit.return_value = False

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
    
    def result_5b(
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
        mock_define_size_of_plot_in_abin,
        mock_show_47b_final_orbit,
        mock_show_instability_limit,
        ):
        result_name = 'result_5b'

        data_ids= [
            data_id('GROUP_sustain_circ_accrete_the_second/low_h_low_acc',  0, 'low h, low acc'),
            data_id('GROUP_sustain_circ_accrete_the_second/low_h_mid_acc',  0, 'low h, mid acc'),

            data_id('GROUP_SC_large_mass_FULL/setup1', 0, r'$40 M_{\oplus}$, setup 2'), #, 'tab:blue'),
            data_id('GROUP_sustain_circ_accrete_the_second/setup1_low_acc', 0, 'setup 1, low acc'),
            data_id('GROUP_sustain_circ_accrete_the_second/setup1_mid_acc', 0, 'setup 1, mid acc'),

            data_id('GROUP_SC_large_mass_FULL/setup2', 0, r'$40 M_{\oplus}$, setup 2'), #, 'tab:blue'),
            data_id('GROUP_sustain_circ_accrete_the_second/setup2_low_acc', 0, 'setup 2, low acc'),
            data_id('GROUP_sustain_circ_accrete_the_second/setup2_mid_acc', 0, 'setup 2, mid acc'),

            data_id('GROUP_sustain_circ_FULL/sustain_circ_setup3_final_ab', 0, 'setup 3, no acc'),
            data_id('GROUP_sustain_circ_accrete/sustain_circ_setup3_final_ab_accrete_0000069', 0, 'setup 3, low acc'),
            data_id('GROUP_sustain_circ_accrete/sustain_circ_setup3_final_ab_accrete_000069', 0, 'setup 3, mid acc'),
            data_id('GROUP_sustain_circ_accrete/sustain_circ_setup3_final_ab_accrete_00069', 0, 'setup 3, large acc'),

            data_id('GROUP_sustain_circ_accrete/sustain_circ_setup4_final_ab_accrete_0000069', 0, 'setup 4, low acc'),
            data_id('GROUP_sustain_circ_accrete/sustain_circ_setup4_final_ab_accrete_000069', 0, 'setup 4, mid acc'),
            data_id('GROUP_sustain_circ_accrete/sustain_circ_setup4_final_ab_accrete_00069', 0, 'setup 4, large acc'),
        ]
        # mock_selectDataToPlot.return_value = data_name
        mock_select_many_data_ids_to_overlay.return_value = data_ids
        # directories = Navigation_helper.Directories(data_name)
        mock_selectDataFileToPlot.return_value = 0
        # mock_selectManyDataFilesToPlot.return_value = [0, 240]
        mock_selectPlottingRange.return_value = [0, 98]
        mock_selectObjectToPlot.return_value = Kep47b_astrophysical_object
        # mock_selectFunctionsToRun.return_value = 'all'
        mock_define_legend_name.return_value = [data_ids[i].legend_name for i in range(len(data_ids))]
        mock_name_the_plot.return_value = ''
        mock_define_save_plot.side_effect = define_new_define_save_plot_fn(result_name)
        mock_select_a_or_e_to_plot.return_value = [Nbody_Characteristics.possible_data_to_plot.semi_major_axis]
        mock_select_object_config_to_plot.return_value = [Kep47b_astrophysical_object] 

        mock_show_47b_final_orbit.return_value = False
        mock_show_instability_limit.return_value = False

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

        # mock_select_averaging_length.return_value = 10
        # Disc_Characteristics.plot_one_multiple_data_sets_gap_parameters_out()

    def result_6a(
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
        mock_define_size_of_plot_in_abin,
        mock_show_47b_final_orbit,
        mock_show_instability_limit,
        ):
        result_name = 'result_6a'
        mock_name_the_plot.return_value = ''

        mock_select_object_config_to_plot.return_value = [Kep47d_astrophysical_object, Kep47b_astrophysical_object] 
        # mock_selectObjectsToPlot.return_value = [2, 3], [Kep47b_astrophysical_object, Kep47d_astrophysical_object] 
        mock_define_size_of_plot_in_abin.return_value = 15.

        mock_show_47b_final_orbit.return_value = False
        mock_show_instability_limit.return_value = False

        mock_select_averaging_length.return_value = 5
        mock_select_a_or_e_to_plot.return_value = [Nbody_Characteristics.possible_data_to_plot.semi_major_axis]

        data = data_id('GROUP_planets_swap/setup1', 0, 'setup1')
        mock_selectPlottingRange.return_value = [0, 141]
        mock_selectDataToPlot.return_value = data.name
        mock_define_save_plot.side_effect = define_new_define_save_plot_fn('{}_{}'.format(result_name, data.legend_name))
        Nbody_Characteristics.plot_one_using_dat_planet_and_cavity()

        data = data_id('GROUP_planets_swap/setup2', 0, 'setup2')
        mock_selectPlottingRange.return_value = [0, 139]
        mock_selectDataToPlot.return_value = data.name
        mock_define_save_plot.side_effect = define_new_define_save_plot_fn('{}_{}'.format(result_name, data.legend_name))
        Nbody_Characteristics.plot_one_using_dat_planet_and_cavity()

    def result_6b(
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
        mock_define_size_of_plot_in_abin,
        mock_show_47b_final_orbit,
        mock_show_instability_limit,
        ):
        result_name = 'result_6b'

        data_ids = [
            data_id('GROUP_can_47d_keep_47b_stable_BETTER_full/setup1',   0, 'setup1'),
            data_id('GROUP_can_47d_keep_47b_stable_BETTER_full/setup2',   0, 'setup2'),
            data_id('GROUP_can_47d_keep_47b_stable_BETTER_full/setup3',   0, 'setup3'),
            data_id('GROUP_can_47d_keep_47b_stable_BETTER_full/setup4',   0, 'setup4'),
        ]
        mock_select_many_data_ids_to_overlay.return_value = data_ids
        # directories = Navigation_helper.Directories(data_name)
        mock_selectDataFileToPlot.return_value = 0
        mock_selectManyDataFilesToPlot.return_value = [210,211,212,213,214,215]
        mock_selectPlottingRange.return_value = [0, 280] # [0, 215]
        # mock_selectFunctionsToRun.return_value = 'all'
        mock_define_legend_name.return_value = [data_ids[i].legend_name for i in range(len(data_ids))]
        mock_name_the_plot.return_value = ''
        mock_define_size_of_plot_in_abin.return_value = 15.

        Kep47b_astrophysical_object.id = 3
        Kep47d_astrophysical_object.id = 2
        
        mock_select_a_or_e_to_plot.return_value = [Nbody_Characteristics.possible_data_to_plot.semi_major_axis]
        mock_select_object_config_to_plot.return_value = [Kep47d_astrophysical_object, Kep47b_astrophysical_object] 
        mock_selectObjectsToPlot.return_value = [2, 3], [Kep47d_astrophysical_object, Kep47b_astrophysical_object] 

        mock_show_47b_final_orbit.return_value = False
        mock_show_instability_limit.return_value = False

        mock_selectObjectToPlot.return_value = Kep47b_astrophysical_object
        mock_define_save_plot.side_effect = define_new_define_save_plot_fn('{}'.format(result_name))


        # mock_select_averaging_length.return_value = 50
        # Nbody_Characteristics.plot_many_data_id_using_dat()

        # mock_select_averaging_length.return_value = 50
        # Disc_Characteristics.plot_one_multiple_data_sets_overlayed_disc_e_avg()   

        # mock_select_averaging_length.return_value = 10
        # Disc_Characteristics.plot_one_multiple_data_sets_gap_parameters_out()


        for dat in data_ids:
            mock_selectDataToPlot.return_value = dat.name
            mock_define_save_plot.side_effect = define_new_define_save_plot_fn('{}_{}'.format(result_name, dat.legend_name))

            mock_select_averaging_length.return_value = 5
            mock_select_a_or_e_to_plot.return_value = [Nbody_Characteristics.possible_data_to_plot.semi_major_axis]
            Nbody_Characteristics.plot_one_using_dat_planet_and_cavity()
            mock_select_a_or_e_to_plot.return_value = [Nbody_Characteristics.possible_data_to_plot.eccentricity]
            Nbody_Characteristics.plot_one_using_dat_planet_and_cavity()
            Nbody_Characteristics.plot_resonance_dat()
            # TwoD_sigma.plot_many()

        
        Kep47b_astrophysical_object.id = 2
        Kep47d_astrophysical_object.id = 3

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
        mock_define_size_of_plot_in_abin,
        mock_show_47b_final_orbit,
        mock_show_instability_limit,
        ):
        result_name = 'result_7'

        # data_ids= [
        #     data_id('GROUP_low_alpha_low_h_FULL/control',           0, r'h=0.03, $\alpha$=$10^{-4}$'),
        #     data_id('GROUP_low_alpha_low_h_FULL/low_alpha',         0, r'h=0.03, $\alpha$=$10^{-5}$'),
        #     data_id('GROUP_low_alpha_low_h_FULL/less_low_h',        0, r'h=0.02, $\alpha$=$10^{-4}$'),
        #     data_id('GROUP_low_alpha_low_h_FULL/low_h',             0, r'h=0.01, $\alpha$=$10^{-4}$'),
        #     data_id('GROUP_low_alpha_low_h_FULL/low_alpha_low_h',   0, r'h=0.01, $\alpha$=$10^{-5}$'),
        # ]

        data_ids= [
            data_id('GROUP_low_alpha_low_h_FULL/control',           0, 'setup 1'),
            data_id('GROUP_low_alpha_low_h_FULL/low_alpha',         0, r'low $\alpha$'),
            data_id('GROUP_low_alpha_low_h_FULL/less_low_h',        0, 'low h'),
            data_id('GROUP_low_alpha_low_h_FULL/low_h',             0, 'very low h'),
            data_id('GROUP_low_alpha_low_h_FULL/low_alpha_low_h',   0, r'very low h, low $\alpha$'),
        ]


        # mock_selectDataToPlot.return_value = data_name
        mock_select_many_data_ids_to_overlay.return_value = data_ids
        # directories = Navigation_helper.Directories(data_name)
        mock_selectDataFileToPlot.return_value = 0
        # mock_selectManyDataFilesToPlot.return_value = [0, 240]
        mock_selectPlottingRange.return_value = [0, 300]
        mock_selectObjectToPlot.return_value = Kep47b_astrophysical_object
        # mock_selectFunctionsToRun.return_value = 'all'
        mock_define_legend_name.return_value = [data_ids[i].legend_name for i in range(len(data_ids))]
        mock_name_the_plot.return_value = ''
        mock_define_save_plot.side_effect = define_new_define_save_plot_fn(result_name)
        mock_select_a_or_e_to_plot.return_value = [Nbody_Characteristics.possible_data_to_plot.semi_major_axis]
        mock_select_object_config_to_plot.return_value = [Kep47b_astrophysical_object] 
        mock_define_size_of_plot_in_abin.return_value = 18.

        mock_show_47b_final_orbit.return_value = True
        mock_show_instability_limit.return_value = False

        mock_select_averaging_length.return_value = 50
        # mock_select_a_or_e_to_plot.return_value = [Nbody_Characteristics.possible_data_to_plot.semi_major_axis]
        # Nbody_Characteristics.plot_many_data_id_using_dat()
        # mock_select_a_or_e_to_plot.return_value = [Nbody_Characteristics.possible_data_to_plot.eccentricity]
        # Nbody_Characteristics.plot_many_data_id_using_dat()
        
        # mock_select_averaging_length.return_value = 50
        # Disc_Characteristics.plot_one_multiple_data_sets_overlayed_disc_e_avg()   

        # mock_select_averaging_length.return_value = 10
        # Disc_Characteristics.plot_one_multiple_data_sets_gap_parameters_out()
        
        low_h_name = 'GROUP_low_alpha_low_h_FULL/low_h'
        mock_selectDataToPlot.return_value = low_h_name
        directories = Navigation_helper.Directories(low_h_name, save_plots_local_to_data=True)
        mock_selectDataFileToPlot.return_value = Navigation_helper.findMaxFileNumber(directories.out_dir)
        mock_define_save_plot.side_effect = define_new_define_save_plot_fn('{}_low_h'.format(result_name))
        TwoD_sigma.plot_one()

        less_low_h_name = 'GROUP_low_alpha_low_h_FULL/less_low_h'
        mock_selectDataToPlot.return_value = less_low_h_name
        directories = Navigation_helper.Directories(less_low_h_name, save_plots_local_to_data=True)
        mock_selectDataFileToPlot.return_value = Navigation_helper.findMaxFileNumber(directories.out_dir)
        mock_define_save_plot.side_effect = define_new_define_save_plot_fn('{}_less_low_h_name'.format(result_name))
        TwoD_sigma.plot_one()
    
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
        mock_define_size_of_plot_in_abin,
        mock_show_47b_final_orbit,
        mock_show_instability_limit,
        ):
        result_name = 'result_8'

        data_ids = [
            data_id('GROUP_late_stage_disc/setup1',   0, 'setup1'),
            data_id('GROUP_late_stage_disc/setup2',   0, 'setup2'),
            data_id('GROUP_late_stage_disc/setup3',   0, 'setup3'),
            data_id('GROUP_late_stage_disc/setup4',   0, 'setup4'),
        ]
        mock_select_many_data_ids_to_overlay.return_value = data_ids
        # directories = Navigation_helper.Directories(data_name)
        mock_selectDataFileToPlot.return_value = 0
        # mock_selectManyDataFilesToPlot.return_value = [0, 240]
        mock_selectPlottingRange.return_value = [0, 10]
        # mock_selectFunctionsToRun.return_value = 'all'
        mock_define_legend_name.return_value = [data_ids[i].legend_name for i in range(len(data_ids))]
        mock_name_the_plot.return_value = ''
        
        mock_select_a_or_e_to_plot.return_value = [Nbody_Characteristics.possible_data_to_plot.semi_major_axis]
        mock_select_object_config_to_plot.return_value = [Kep47b_astrophysical_object,] 
        mock_selectObjectsToPlot.return_value = [2], [Kep47b_astrophysical_object] 
        mock_define_size_of_plot_in_abin.return_value = 15.

        mock_show_47b_final_orbit.return_value = False
        mock_show_instability_limit.return_value = False

        mock_define_save_plot.side_effect = define_new_define_save_plot_fn('{}'.format(result_name))
        mock_selectObjectToPlot.return_value = Kep47b_astrophysical_object

        mock_select_averaging_length.return_value = 50
        Nbody_Characteristics.plot_many_data_id_using_dat()
        

if __name__ == '__main__':
    # report_plots.intro_disc_distribution()
    # report_plots.theory()

    report_plots.result_1_LUKE()
    report_plots.result_2_LUKE()
    report_plots.result_3_LUKE()

    # report_plots.result_3()
    # report_plots.result_4b()
    # report_plots.result_5b()
    # report_plots.result_6a()
    # report_plots.result_6b()
    report_plots.result_7()
    # report_plots.result_8()

    # generate_all_plots()