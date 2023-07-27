import math
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Union



from dynamotable.convention import COLUMN_NAMES
from dynamotable.table_map import table_map_read
from dynamotable.dynamotable import read
from dynamotable.dynamotable import write

def add_columns_for_duplicate_exclusion(table_dataframe):
    table_dataframe['absolute_x'] = table_dataframe['dx'].to_numpy() + table_dataframe['x'].to_numpy()
    table_dataframe['absolute_y'] = table_dataframe['dy'].to_numpy() + table_dataframe['y'].to_numpy()
    table_dataframe['absolute_z'] = table_dataframe['dz'].to_numpy() + table_dataframe['z'].to_numpy()
    table_dataframe['points_on_lattice'] = 0
    table_dataframe['to_delete'] = 0
    return table_dataframe


def split_dataframe_according_to_mt_number(table_dataframe):
    grouped_table_dataframe = table_dataframe.groupby(table_dataframe['annotation'])
    list_of_table_dataframe = []

    for i in range(1,(table_dataframe['annotation'].to_numpy().max()) + 1):
        list_of_table_dataframe.append(grouped_table_dataframe.get_group(i))
    return list_of_table_dataframe


def distance_within_range(array2,
                          array1,
                          limit):
    distance = math.sqrt((array2[0] - array1[0]) ** 2 + (array2[1] - array1[1]) ** 2 + (array2[2] - array1[2]) ** 2)
    if distance <= limit and array2[3] > array1[3]:
        return 1
    elif distance <= limit and array2[3] == array1[3] and array2[4] <= array1[4]:
        return 1
    else:
        return 0


def neighbor_with_better_geometry(array_of_point1,
                                  array_of_other_points,
                                  limit):
    array_of_neighbor_count = np.apply_along_axis(distance_within_range,
                                                  axis=1,
                                                  arr=array_of_other_points,
                                                  array1=array_of_point1,
                                                  limit=limit)
    if np.sum(array_of_neighbor_count) > 1:
        return 1
    else:
        return 0


def neighbor_with_better_geometry_in_table_dataframe(table_dataframe,
                                                     limit):
    array_for_neighbor_distance = table_dataframe[['absolute_x', 'absolute_y', 'absolute_z', 'points_on_lattice', 'tag']].to_numpy()
    array_of_neighbor_count = np.apply_along_axis(neighbor_with_better_geometry,
                                                  axis=1,
                                                  arr=array_for_neighbor_distance,
                                                  array_of_other_points=array_for_neighbor_distance,
                                                  limit=limit)
    table_dataframe.loc[:,'to_delete'] = array_of_neighbor_count
    return table_dataframe


def distance_within_two_ranges(array2,
                               array1,
                               lower_limit_first_bin, upper_limit_first_bin,
                               lower_limit_second_bin, upper_limit_second_bin):
    distance = math.sqrt((array2[0] - array1[0]) ** 2 + (array2[1] - array1[1]) ** 2 + (array2[2] - array1[2]) ** 2)
    in_distance = lower_limit_first_bin <= distance <= upper_limit_first_bin or lower_limit_second_bin <= distance <= upper_limit_second_bin
    if in_distance:
        return 1
    else:
        return 0

def number_of_neighbors_within_two_ranges(array_of_point1,
                                          array_of_other_points,
                                          lower_limit_first_bin, upper_limit_first_bin,
                                          lower_limit_second_bin, upper_limit_second_bin):
    array_of_neighbor_count = np.apply_along_axis(distance_within_two_ranges,
                                                  axis=1,
                                                  arr=array_of_other_points,
                                                  array1=array_of_point1,
                                                  lower_limit_first_bin=lower_limit_first_bin,
                                                  upper_limit_first_bin=upper_limit_first_bin,
                                                  lower_limit_second_bin=lower_limit_second_bin,
                                                  upper_limit_second_bin=upper_limit_second_bin)
    return np.sum(array_of_neighbor_count)

def number_of_neighbors_in_table_dataframe_within_two_ranges(table_dataframe,
                                                             lower_limit_first_bin, upper_limit_first_bin,
                                                             lower_limit_second_bin, upper_limit_second_bin):
    array_for_neighbor_distance = table_dataframe[['absolute_x', 'absolute_y', 'absolute_z']].to_numpy()
    array_of_neighbor_count = np.apply_along_axis(number_of_neighbors_within_two_ranges,
                                                  axis=1,
                                                  arr=array_for_neighbor_distance,
                                                  array_of_other_points=array_for_neighbor_distance,
                                                  lower_limit_first_bin=lower_limit_first_bin,
                                                  upper_limit_first_bin=upper_limit_first_bin,
                                                  lower_limit_second_bin=lower_limit_second_bin,
                                                  upper_limit_second_bin=upper_limit_second_bin)
    table_dataframe.loc[:,'points_on_lattice'] = array_of_neighbor_count
    return table_dataframe


def remove_identified_duplicates_and_intermediate_columns(table_dataframe):
    to_delete = table_dataframe[table_dataframe['to_delete'] == 1].index
    table_dataframe.drop(index=to_delete, inplace=True)
    table_dataframe.drop(columns=['absolute_x', 'absolute_y', 'absolute_z', 'points_on_lattice', 'to_delete'], inplace=True)
    return table_dataframe


def remove_duplicates_full(table_name,
                           lower_limit_first_bin, upper_limit_first_bin,
                           lower_limit_second_bin, upper_limit_second_bin,
                           proximity_limit):
    list_of_dataframes = split_dataframe_according_to_mt_number(add_columns_for_duplicate_exclusion(read(table_name)))
    list_of_updated_dataframes = []
    for dataframe in list_of_dataframes:
        list_of_updated_dataframes.append(
            remove_identified_duplicates_and_intermediate_columns(
                neighbor_with_better_geometry_in_table_dataframe(
                    number_of_neighbors_in_table_dataframe_within_two_ranges(
                        dataframe,
                        lower_limit_first_bin, upper_limit_first_bin,
                        lower_limit_second_bin, upper_limit_second_bin),
                    proximity_limit)))
    return pd.concat(list_of_updated_dataframes)


#print(remove_duplicates_full('test_table.tbl', 1.5, 2.5, 3.5, 4.5, 1.5))
write(remove_duplicates_full('converted_table.tbl', 4.5, 5.05, 9.0, 9.5, 2.2), "output.tbl")
#input_table_file = 'test_table.tbl'
#table = read(input_table_file)
#table = add_columns_for_duplicate_exclusion(read(input_table_file))
#table_list = split_dataframe_according_to_mt_number(table)
#first_table = table_list[0]
#second_table = table_list[1]
#third_table = table_list[2]

#fourth_table = number_of_neighbors_in_table_dataframe_within_two_ranges(third_table, 1.5, 2.5, 3.5, 4.5)
#fifth_table = neighbor_with_better_geometry_in_table_dataframe(fourth_table, 1.5)
#sixth_table = remove_identified_duplicates_and_intermediate_columns(fifth_table)
#print(fourth_table)
#print(fifth_table)
#print(sixth_table)

