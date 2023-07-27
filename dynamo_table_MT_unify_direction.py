import math
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Union


from dynamotable.convention import COLUMN_NAMES
from dynamotable.table_map import table_map_read
from dynamotable.dynamotable import read
from dynamotable.dynamotable import write

def split_dataframe_according_to_mt_number(table_dataframe):
    grouped_table_dataframe = table_dataframe.groupby(table_dataframe['annotation'])
    list_of_table_dataframe = []
    for i in range(1,(table_dataframe['annotation'].to_numpy().max()) + 1):
        list_of_table_dataframe.append(grouped_table_dataframe.get_group(i))
    return list_of_table_dataframe


def test_for_tilt_difference(table_dataframe_orig, table_dataframe_aligned):
    tdrot_orig = table_dataframe_orig['tdrot'].to_numpy()
    tdrot_aligned = table_dataframe_aligned['tdrot'].to_numpy()
    accumulated_angle_difference = 0
    for i in range(0,len(tdrot_orig)):
        if tdrot_orig[i] < 0:
            tdrot_orig[i] = tdrot_orig[i] + 360
        if tdrot_aligned[i] < 0:
            tdrot_aligned[i] = tdrot_aligned[i] + 360
        if abs(tdrot_orig[i] - tdrot_aligned[i]) > 90:
            accumulated_angle_difference += 1
    return (accumulated_angle_difference / len(tdrot_orig))



def unify_direction_per_filament(table_dataframe_orig, table_dataframe_aligned, table_data_frame_to_be_updated):
    average_angle_difference = test_for_tilt_difference(table_dataframe_orig, table_dataframe_aligned)
    if 0.33 < average_angle_difference < 0.66:
        print(f"check filament #{table_dataframe_orig.iloc[1]['annotation']}, {round(average_angle_difference,2)} of particles support flipping")
    if 0.5 < average_angle_difference:
        table_data_frame_to_be_updated.loc[:, 'tilt'] = round(table_data_frame_to_be_updated.loc[:, 'tilt'] + 180, 2)
        print(f"filament #{table_dataframe_orig.iloc[1]['annotation']} was flipped")
        inv_table_data_frame_to_be_updated = table_data_frame_to_be_updated.loc[::-1] #reverses Order so that all MT start on the same end
        table_data_frame_to_be_updated.loc[:, 'tag'] = inv_table_data_frame_to_be_updated.loc[:, 'tag'] #then updates particle numbers
    return table_data_frame_to_be_updated


def unify_direction(table_dataframe_orig, table_dataframe_aligned, table_data_frame_to_be_updated):
    list_table_data_frame_orig = split_dataframe_according_to_mt_number(table_dataframe_orig)
    list_table_data_frame_aligned = split_dataframe_according_to_mt_number(table_dataframe_aligned)
    list_table_data_frame_to_be_updated = split_dataframe_according_to_mt_number(table_data_frame_to_be_updated)
    list_table_data_frame_unified = []
    for i in range(0,len(list_table_data_frame_orig)):
        list_table_data_frame_unified.append(
            unify_direction_per_filament(list_table_data_frame_orig[i], list_table_data_frame_aligned[i], list_table_data_frame_to_be_updated[i]))
    return pd.concat(list_table_data_frame_unified, ignore_index=True) #Resets indices to preserve order just in case this will be used later


#print(remove_duplicates_full('test_table.tbl', 1.5, 2.5, 3.5, 4.5, 1.5))
#write(remove_duplicates_full('converted_table.tbl', 4.5, 5.05, 9.0, 9.5, 2.2), "output.tbl")
input_table_file1 = 'finding_position_ref_run001_ite1.tbl'
input_table_file3 = 'finding_position_ref_run002_ite1.tbl'
input_table_file5 = 'initial_positions_bin4_crop.tbl'
table1 = read(input_table_file1)
table3 = read(input_table_file3)
table5 = read(input_table_file5)
write(unify_direction(table1, table3, table5),'finding_position_ref_unified_direction.tbl')
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

