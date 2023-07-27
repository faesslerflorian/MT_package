# Import pandas package
import pandas as pd
import numpy as np
import math
import argparse

from scipy.interpolate import CubicSpline
from scipy.spatial.transform import Rotation
from dynamotable_autofill import read
from dynamotable.dynamotable import write
from dynamotable.convention import COLUMN_NAMES


def split_dataframe_according_to_mt_number(table_dataframe):
    grouped_table_dataframe = table_dataframe.groupby(table_dataframe['annotation'])
    list_of_table_dataframe = []
    for i in range(1,(table_dataframe['annotation'].to_numpy().max()) + 1):
        list_of_table_dataframe.append(grouped_table_dataframe.get_group(i).reset_index(drop=True))
    return list_of_table_dataframe


def geo_vector(point1,
               point2):
    x = point2[0] - point1[0]
    y = point2[1] - point1[1]
    z = point2[2] - point1[2]
    v = (x, y, z)
    return v


def table_data_frame_to_spline_interpolated_contour(table_dataframe,
                                                    step):
    x_values = table_dataframe['dx'].to_numpy() + table_dataframe['x'].to_numpy()
    y_values = table_dataframe['dy'].to_numpy() + table_dataframe['y'].to_numpy()
    z_values = table_dataframe['dz'].to_numpy() + table_dataframe['z'].to_numpy()
    t_values = np.array([])
    current_t_value = 0
    t_values = np.append(t_values, current_t_value)
    interpolated_contour = []
    for i in range(1, len(x_values)):
        current_t_value = current_t_value + math.sqrt((x_values[i] - x_values[i - 1]) ** 2
                                                      + (y_values[i] - y_values[i - 1]) ** 2
                                                      + (z_values[i] - z_values[i - 1]) ** 2)
        t_values = np.append(t_values, current_t_value)
    cx = CubicSpline(t_values, x_values)
    cy = CubicSpline(t_values, y_values)
    cz = CubicSpline(t_values, z_values)
    new_t_values = np.arange(0, math.floor(t_values[-1]), step)
    new_x_values = cx(new_t_values)
    new_y_values = cy(new_t_values)
    new_z_values = cz(new_t_values)
    for idx, point in enumerate(new_t_values):
        interpolated_contour.append([round(new_x_values[idx], 2),
                                     round(new_y_values[idx], 2),
                                     round(new_z_values[idx], 2)])
    return interpolated_contour


def length_three_d_vector(vector):
    length = math.sqrt(vector[0] ** 2 + vector[1] ** 2 + vector[2] ** 2)
    return length


def normalize_vector(vector):
    length = length_three_d_vector(vector)
    normalized_vector = (vector[0]/length, vector[1]/length, vector[2]/length)
    return normalized_vector


def signed_angle_between_orientation_vectors(vector1,
                                             vector2,
                                             normal):
    normalized_normal = normalize_vector(normal)
    signed_angle = math.degrees(math.atan2(np.dot(np.cross(vector1, vector2), normalized_normal),
                                           np.dot(vector1, vector2)))
    return signed_angle


def rotate_vector_around_axis(initial_vector,
                              rotation_axis,
                              rotation_radians):
    rotation_vector = [rotation_radians * rotation_axis[0], rotation_radians * rotation_axis[1], rotation_radians * rotation_axis[2]]
    rotation = Rotation.from_rotvec(rotation_vector)
    rotated_vector = rotation.apply(initial_vector)
    return rotated_vector


def tdrot_and_tilt_for_z_alignment(vector):
    tdrot = -round(math.degrees(math.atan2(-vector[0], vector[1])),2)
    if vector[0] == vector[1] == 0:
        if vector[2] > 0:
            tilt = 0
        else:tilt = 180
    else:
        normalized_vector = normalize_vector(vector)
        tilt = math.degrees(math.acos(normalized_vector[2]))
    return [tdrot, tilt]


def apply_tdrot_and_tilt_to_x(tdrot,
                              tilt):
    rotation_to_apply = Rotation.from_euler('ZX',[tdrot, tilt], degrees=True)
    rotated_x = rotation_to_apply.apply([1, 0, 0])
    return [rotated_x[0], -rotated_x[1], rotated_x[2]]


def all_eulers_for_z_and_x_alignment(direction_vector_sampling_point,
                                     direction_vector_axis):
    [tdrot, tilt] = tdrot_and_tilt_for_z_alignment(direction_vector_sampling_point)
    intermediate_x_vector = apply_tdrot_and_tilt_to_x(tdrot, tilt)
    narot = signed_angle_between_orientation_vectors(direction_vector_axis,
                                                     intermediate_x_vector,
                                                     direction_vector_sampling_point)
    return [round(tdrot, 2), round(tilt, 2), round(narot, 2)]


def sample_tube_at_one_position(point,
                                direction_vector,
                                tube_radius,
                                sampling_rate_circle):
    sampling_angle = 360 / sampling_rate_circle
    sampling_radians = math.radians(sampling_angle)
    normalized_direction_vector = normalize_vector(direction_vector)
    [tdrot_current_vector, tilt_current_vector] = tdrot_and_tilt_for_z_alignment(direction_vector)
    new_z_after_first_2_eulers = apply_tdrot_and_tilt_to_x(tdrot_current_vector, tilt_current_vector)
    direction_vector_current_sampling_point = [new_z_after_first_2_eulers[0] * tube_radius,
                                               new_z_after_first_2_eulers[1] * tube_radius,
                                               new_z_after_first_2_eulers[2] * tube_radius]
    eulers_current_vector = all_eulers_for_z_and_x_alignment(direction_vector_current_sampling_point, direction_vector)
    coordinates_and_direction_vectors_of_all_sampling_points = [[[round(point[0] + direction_vector_current_sampling_point[0],2),
                                                                  round(point[1] + direction_vector_current_sampling_point[1],2),
                                                                  round(point[2] + direction_vector_current_sampling_point[2],2)],
                                                                 eulers_current_vector]]
    for i in range(1, sampling_rate_circle, 1):
        direction_vector_current_sampling_point = rotate_vector_around_axis(direction_vector_current_sampling_point,
                                                                            normalized_direction_vector,
                                                                            sampling_radians)
        eulers_current_vector = all_eulers_for_z_and_x_alignment(direction_vector_current_sampling_point, direction_vector)
        coordinates_and_direction_vectors_of_all_sampling_points.append([[round(point[0] + direction_vector_current_sampling_point[0], 2),
                                                                          round(point[1] + direction_vector_current_sampling_point[1], 2),
                                                                          round(point[2] + direction_vector_current_sampling_point[2], 2)],
                                                                         eulers_current_vector])
    return coordinates_and_direction_vectors_of_all_sampling_points


def sample_tube(table_dataframe,
                sampling_rate_contour,
                tube_radius,
                sampling_rate_circle):
    resampled_contour = table_data_frame_to_spline_interpolated_contour(table_dataframe, sampling_rate_contour)
    sampling_points_on_tube = []
    for point in range(0, len(resampled_contour)):
        if point == 0:
            sampling_points_on_tube += sample_tube_at_one_position(resampled_contour[0],
                                                                   geo_vector(resampled_contour[0],
                                                                              resampled_contour[1]),
                                                                   tube_radius,
                                                                   sampling_rate_circle)
        else:
            sampling_points_on_tube += sample_tube_at_one_position(resampled_contour[point],
                                                                   geo_vector(resampled_contour[point - 1],
                                                                              resampled_contour[point]),
                                                                   tube_radius, sampling_rate_circle)
    return sampling_points_on_tube


def tube_surface_dataframe_from_tube_center_dataframe(table_dataframe,
                                                      sampling_rate_contour,
                                                      tube_radius,
                                                      sampling_rate_circle):
    sampling_points_on_tube = sample_tube(table_dataframe,
                                          sampling_rate_contour,
                                          tube_radius,
                                          sampling_rate_circle)
    tube_surface_dataframe = pd.DataFrame([], columns=COLUMN_NAMES)
    number_of_sampling_points = len(sampling_points_on_tube)
    number_of_table_columns = len(COLUMN_NAMES)
    common_data_on_tube = table_dataframe.iloc[0].to_numpy(copy=True)
    tube_surface_array = np.full((number_of_sampling_points,number_of_table_columns),common_data_on_tube)
    tube_surface_array[:, 1] = 1
    tube_surface_array[:, [2, 3, 4, 5, 9, 10, 11, 17, 18, 26, 27, 28, 29, 30, 31, 33, 34, 40, 41]] = 0
    x_values = np.empty(number_of_sampling_points)
    y_values = np.empty(number_of_sampling_points)
    z_values = np.empty(number_of_sampling_points)
    tdrot_values = np.empty(number_of_sampling_points)
    tilt_values = np.empty(number_of_sampling_points)
    narot_values = np.empty(number_of_sampling_points)
    for idx, sampling_point in enumerate(sampling_points_on_tube):
        coords = sampling_point[0]
        x_values[idx] = coords[0]
        y_values[idx] = coords[1]
        z_values[idx] = coords[2]
        eulers = sampling_point[1]
        tdrot_values[idx] = eulers[0]
        tilt_values[idx] = eulers[1]
        narot_values[idx] = eulers[2]
    tube_surface_array[:, 6] = tdrot_values
    tube_surface_array[:, 7] = tilt_values
    tube_surface_array[:, 8] = narot_values
    tube_surface_array[:, 23] = x_values
    tube_surface_array[:, 24] = y_values
    tube_surface_array[:, 25] = z_values
    tube_surface_dataframe = pd.DataFrame(tube_surface_array, columns=COLUMN_NAMES)
    return tube_surface_dataframe


def complete_surface_dataframe_from_complete_center_dataframe(complete_center_dataframe,
                                                              sampling_rate_contour,
                                                              tube_radius,
                                                              sampling_rate_circle):
    list_of_tube_center_dataframes = split_dataframe_according_to_mt_number(complete_center_dataframe)
    list_of_tube_surface_dataframes = []
    for tube_center_dataframe in list_of_tube_center_dataframes:
        list_of_tube_surface_dataframes.append(tube_surface_dataframe_from_tube_center_dataframe(tube_center_dataframe,
                                                                                                 sampling_rate_contour,
                                                                                                 tube_radius,
                                                                                                 sampling_rate_circle))
    complete_surface_dataframe = pd.concat(list_of_tube_surface_dataframes,
                                           axis=0,
                                           ignore_index=True)
    complete_surface_dataframe.loc[:, 'tag'] = np.arange(1, complete_surface_dataframe.shape[0] + 1)
    return complete_surface_dataframe


def complete_surface_table_from_complete_center_table(complete_center_table,
                                                      complete_surface_table,
                                                      sampling_rate_contour,
                                                      tube_radius,
                                                      sampling_rate_circle):
    complete_center_dataframe = read(complete_center_table)
    complete_surface_dataframe = complete_surface_dataframe_from_complete_center_dataframe(complete_center_dataframe,
                                                                                           sampling_rate_contour,
                                                                                           tube_radius,
                                                                                           sampling_rate_circle)
    write(complete_surface_dataframe, complete_surface_table)
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='samples surfaces of tubes whose center lines are represented in '
                                                 'a dynamo .tbl, tubes are identified by column 23/"annotation"'
                                                 'values have to be provided in the same scale as in the .tbl'
                                                 'output is a dynamo .tbl as well')
    parser.add_argument('input',
                        type=str,
                        help='dynamo .tbl used as input')
    parser.add_argument('--tube_length_sampling_point_distance',
                        default=1,
                        type=float,
                        help='defines defines the distance between sampling points along the filament/tube axis'
                             'default = 1')
    parser.add_argument('--tube_radius',
                        default=0.1,
                        type=float,
                        help='defines in which radius around the center line sampling will be performed.'
                             'default 0.1')
    parser.add_argument('--tube_circumference_sampling_points',
                        default=1,
                        type=int,
                        help='defines how many regularly spaced points on an assumed tube circumference'
                             'are sampled. '
                             'default = 1')
    parser.add_argument('--output',
                        default='output.tbl',
                        type=str,
                        help='defines the name of the .tbl that is written out.'
                             'default = output.tbl')

    args = parser.parse_args()

    complete_surface_table_from_complete_center_table(args.input,
                                                      args.output,
                                                      args.tube_length_sampling_point_distance,
                                                      args.tube_radius,
                                                      args.tube_circumference_sampling_points)
