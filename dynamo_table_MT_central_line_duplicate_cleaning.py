import pandas as pd
import numpy as np
import math
import argparse

from scipy.interpolate import CubicSpline
from scipy.spatial.transform import Rotation
from dynamotable_autofill import read
from dynamotable.dynamotable import write
from dynamotable.convention import COLUMN_NAMES
from dynamo_table_dataframe_functions import split_table_dataframe_according_to_mt_number
from dynamo_table_dataframe_functions import array_for_geometry_from_table_dataframe
from dynamo_table_dataframe_functions import array_for_distance_from_table_dataframe
from dynamo_table_dataframe_functions import apply_euler_angle_rotations_to_vector


# Assumes alignment of the z-axis with the central line of the MT and the x-axis with the "radius" at the current position
def project_central_line_position_to_lattice(one_d_array_for_geometry, tube_radius):
    [tag, tdrot, tilt, narot, absolut_x, absolut_y, absolut_z] = one_d_array_for_geometry
    shift_vector = apply_euler_angle_rotations_to_vector(tdrot, tilt, narot, [tube_radius, 0, 0])
    updated_absolut_x = round((absolut_x + shift_vector[0]), 2)
    updated_absolut_y = round((absolut_y + shift_vector[1]), 2)
    updated_absolut_z = round((absolut_z + shift_vector[2]), 2)
    updated_tdrot = round(tdrot, 2)
    updated_tilt = round((tilt + 90), 2)
    updated_narot = round((narot + 90), 2)
    return np.array([tag, updated_tdrot, updated_tilt, updated_narot, updated_absolut_x, updated_absolut_y, updated_absolut_z])


def project_central_line_positions_to_lattice(two_d_array_for_geometry, tube_radius):
    print(two_d_array_for_geometry)
    output_two_d_array_for_geometry = np.empty(shape=7)
    for one_d_array_for_geometry in two_d_array_for_geometry:
        output_one_d_array_for_geometry = project_central_line_position_to_lattice(one_d_array_for_geometry, tube_radius)
        output_two_d_array_for_geometry = np.vstack((output_two_d_array_for_geometry, output_one_d_array_for_geometry))
    return np.delete(output_two_d_array_for_geometry, 0, 0)


def dynamo_table_dataframe_project_central_line_positions_to_lattice(table_dataframe, tube_radius):
    input_two_d_array_for_geometry = array_for_geometry_from_table_dataframe(table_dataframe)
    output_two_d_array_for_geometry = project_central_line_positions_to_lattice(input_two_d_array_for_geometry, tube_radius)
    table_dataframe.loc[:, 'tag'] = output_two_d_array_for_geometry[:, 0]
    table_dataframe.loc[:, 'tdrot'] = output_two_d_array_for_geometry[:, 1]
    table_dataframe.loc[:, 'tilt'] = output_two_d_array_for_geometry[:, 2]
    table_dataframe.loc[:, 'narot'] = output_two_d_array_for_geometry[:, 3]
    table_dataframe.loc[:, 'x'] = output_two_d_array_for_geometry[:, 4]
    table_dataframe.loc[:, 'y'] = output_two_d_array_for_geometry[:, 5]
    table_dataframe.loc[:, 'z'] = output_two_d_array_for_geometry[:, 6]
    table_dataframe.loc[:, 'dx'] = 0
    table_dataframe.loc[:, 'dy'] = 0
    table_dataframe.loc[:, 'dz'] = 0

    return table_dataframe


def test_2_points_for_proximity(one_d_array_for_distance_1, one_d_array_for_distance_2, distance_cut_off):
    [tag_1, absolut_x_1, absolut_y_1, absolut_z_1] = one_d_array_for_distance_1
    [tag_2, absolut_x_2, absolut_y_2, absolut_z_2] = one_d_array_for_distance_2
    return distance_cut_off >= math.sqrt((absolut_x_1-absolut_x_2) ** 2
                                         + (absolut_y_1-absolut_y_2) ** 2
                                         + (absolut_z_1-absolut_z_2) ** 2)


def test_1_point_and_mt_for_proximity(one_d_array_for_distance_1, two_d_array_for_distance, distance_cut_off):
    for one_d_array_for_distance_2 in two_d_array_for_distance:
        if test_2_points_for_proximity(one_d_array_for_distance_1, one_d_array_for_distance_2, distance_cut_off):
            return True
    return False


def distance_clean_array_for_distance(two_d_array_for_distance, distance_cut_off):
    output_two_d_array_for_distance = np.array([[0,-1000000,-1000000, -1000000]])
    #This is hardcoded to be outside any expected coordinate system! Adjust here in case of missing particles
    for one_d_array_for_distance in two_d_array_for_distance:
        if not test_1_point_and_mt_for_proximity(one_d_array_for_distance, output_two_d_array_for_distance, distance_cut_off):
            output_two_d_array_for_distance = np.vstack((one_d_array_for_distance, output_two_d_array_for_distance))
    output_two_d_array_for_distance = np.flipud(output_two_d_array_for_distance)
    return np.delete(output_two_d_array_for_distance, 0, 0)


def distance_clean_single_mt_table_dataframe(table_dataframe, distance_cut_off):
    two_d_array_for_remaining_particles = distance_clean_array_for_distance(
        array_for_distance_from_table_dataframe(table_dataframe),
        distance_cut_off)
    remaining_tags = two_d_array_for_remaining_particles[:, 0]
    return table_dataframe[table_dataframe['tag'].isin(remaining_tags)]


def distance_clean_complete_table_dataframe(table_dataframe, distance_cut_off):
    single_mt_table_dataframe_list = split_table_dataframe_according_to_mt_number(table_dataframe)
    cleaned_single_mt_table_dataframe_list = []
    for single_MT_table_dataframe in single_mt_table_dataframe_list:
        cleaned_single_mt_table_dataframe_list.append(distance_clean_single_mt_table_dataframe(single_MT_table_dataframe, distance_cut_off))
    return pd.concat(cleaned_single_mt_table_dataframe_list, axis=0, ignore_index=True)




#a = read('column_28_set_to_0.tbl')
#b = dynamo_table_dataframe_project_central_line_positions_to_lattice(a, 12)
#write(b, 'column_28_set_to_0_shifted.tbl')
#c = distance_clean_complete_table_dataframe(b, 2)
#write(c, 'column_28_set_to_0_shifted_cleaned.tbl')



