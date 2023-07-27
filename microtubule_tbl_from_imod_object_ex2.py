import numpy as np
import math
from scipy.interpolate import CubicSpline
from scipy.spatial.transform import Rotation
import argparse
from sympy import Plane, Point3D

def geo_vector(point1, point2):
    x = point2[0] - point1[0]
    y = point2[1] - point1[1]
    z = point2[2] - point1[2]
    v = (x, y, z)
    return v


def length_three_d_vector(vector):
    length = math.sqrt(vector[0] ** 2 + vector[1] ** 2 + vector[2] ** 2)
    return length


def normalize_vector(vector):
    length = length_three_d_vector(vector)
    normalized_vector = (vector[0]/length, vector[1]/length, vector[2]/length)
    return normalized_vector


def signed_angle_between_orientation_vectors(vector1, vector2, normal):
    normalized_normal = normalize_vector(normal)
    signed_angle = math.degrees(math.atan2(np.dot(np.cross(vector1, vector2),normalized_normal), np.dot(vector1, vector2)))
    return signed_angle


def angle_between_orientation_vectors(vector1, vector2):
    an = math.sqrt(vector1[0] * vector1[0] + vector1[1] * vector1[1] + vector1[2] * vector1[2])
    bn = math.sqrt(vector2[0] * vector2[0] + vector2[1] * vector2[1] + vector2[2] * vector2[2])
    ax, ay, az = vector1[0] / an, vector1[1] / an, vector1[2] / an
    bx, by, bz = vector2[0] / bn, vector2[1] / bn, vector2[2] / bn
    return math.degrees(math.acos(ax*bx + ay*by + az*bz))

# reads a text file containing a model and creates a list of entries. For each point, for entries are created:
# first a number to identify the according contour number and 3 numbers to identy x-, y- and z-coordinates


def file_to_list(input_file):
    output_list = []
    for line in open(input_file):
        output_list.extend(line.split())
    return output_list


# uses the list produced by file_to_list to build a list of named points, which are represented by quartuple
# consisting of a contour identifier number and coordinates in x, y and z


def element_list_to_named_point_list(input_list):
    building_named_point = []
    building_named_point_list = []
    counter = 0
    for ele in input_list:
        if counter == 0:
            building_named_point.append(float(ele))
            counter = (counter + 1)
        elif 0 < counter < 3:
            building_named_point.append(float(ele))
            counter = (counter + 1)
        else:
            building_named_point.append(float(ele))
            building_named_point_list.append(building_named_point)
            counter = 0
            building_named_point = []
    return building_named_point_list


# uses the list produced by element_list_to_named_point_list to build a model
# consisting of a list of contours identifier number and coordinates in x, y and z
# consisting of a list of points
# consisting of triples of x, y and z coordinates


def named_point_list_to_model(input_list):
    building_contour = []
    building_model = []
    current_contour = input_list[0]
    current_contour_number = current_contour[0]
    for ele_index, ele in enumerate(input_list):
        if ele[0] == current_contour_number:
            building_contour.append((ele[1], ele[2], ele[3]))
        else:
            building_model.append(building_contour)
            building_contour = [(ele[1], ele[2], ele[3])]
            if ele_index != len(input_list) - 1:
                current_contour = input_list[(ele_index)]
                current_contour_number = current_contour[0]
    building_model.append(building_contour)
    return building_model


# reads a text file containing a model and builds a model in python
# a list of contours, which are a list of points, which are represented by triples)
def file_to_model(file_name):
    return named_point_list_to_model(element_list_to_named_point_list(file_to_list(file_name)))


#Aplly a scaling factor to all points in a contour
def scale_contour(contour, scale):
    output_contour = []
    for point in contour:
        new_point = [0, 0, 0]
        new_point[0] = point[0] * scale
        new_point[1] = point[1] * scale
        new_point[2] = point[2] * scale
        output_contour.append(new_point)
    return output_contour


#Interpolates contours and enforces the respective spacing between the points of the contour
def contour_to_spline_interpolated_contour(contour, step):
    x_values = np.array([])
    y_values = np.array([])
    z_values = np.array([])
    t_values = np.array([])
    current_t_value = 0
    interpolated_contour = []
    for idx, point in enumerate(contour):
        current_x_value = point[0]
        current_y_value = point[1]
        current_z_value = point[2]
        x_values = np.append(x_values, current_x_value)
        y_values = np.append(y_values, current_y_value)
        z_values = np.append(z_values, current_z_value)
        if idx > 0:
            previous_point = contour[idx-1]
            previous_x_value = previous_point[0]
            previous_y_value = previous_point[1]
            previous_z_value = previous_point[2]
            current_t_value = current_t_value + math.sqrt((current_x_value - previous_x_value) ** 2
                                                          + (current_y_value - previous_y_value) ** 2
                                                          + (current_z_value - previous_z_value) ** 2)
        t_values = np.append(t_values, current_t_value)
    cx = CubicSpline(t_values, x_values)
    cy = CubicSpline(t_values, y_values)
    cz = CubicSpline(t_values, z_values)
    new_t_values = np.arange(0, math.floor(t_values[-1]), step)
    new_x_values = cx(new_t_values)
    new_y_values = cy(new_t_values)
    new_z_values = cz(new_t_values)
    for idx, point in enumerate(new_t_values):
        interpolated_contour.append([round(new_x_values[idx],2), round(new_y_values[idx],2), round(new_z_values[idx],2)])
    return interpolated_contour


#Applies contour_to_spline_interpolated_contour to all contours of a model
def model_to_spline_interpolated_model(model, step):
    interpolated_model = []
    for contour in model:
        current_contour = contour_to_spline_interpolated_contour(contour, step)
        interpolated_model.append(current_contour)
    return interpolated_model


def rotate_vector_around_axis(initial_vector, rotation_axis, rotation_radians):
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


def apply_tdrot_and_tilt_to_x(tdrot, tilt):
    rotation_to_apply = Rotation.from_euler('ZX',[tdrot, tilt], degrees=True)
    rotated_x = rotation_to_apply.apply([1, 0, 0])
    return [rotated_x[0], -rotated_x[1], rotated_x[2]]


def all_eulers_for_z_and_x_alignment(direction_vector_sampling_point, direction_vector_axis):
    [tdrot, tilt] = tdrot_and_tilt_for_z_alignment(direction_vector_sampling_point)
    intermediate_x_vector = apply_tdrot_and_tilt_to_x(tdrot, tilt)
    narot = signed_angle_between_orientation_vectors(direction_vector_axis, intermediate_x_vector, direction_vector_sampling_point)
    #narot = angle_between_orientation_vectors(normalize_vector(intermediate_x_vector), normalize_vector(direction_vector_axis))
    return [round(tdrot,2), round(tilt,2), round(narot,2)]


def sample_tube_at_one_position(point, direction_vector, tube_radius, sampling_rate_circle):
    sampling_angle = 360 / sampling_rate_circle
    sampling_radians = math.radians(sampling_angle)
    normalized_direction_vector = normalize_vector(direction_vector)
    [tdrot_current_vector, tilt_current_vector] = tdrot_and_tilt_for_z_alignment(direction_vector)
    new_z_after_first_2_eulers = apply_tdrot_and_tilt_to_x(tdrot_current_vector, tilt_current_vector)
    direction_vector_current_sampling_point = [new_z_after_first_2_eulers[0] * tube_radius, new_z_after_first_2_eulers[1] * tube_radius, new_z_after_first_2_eulers[2]* tube_radius]
    eulers_current_vector = all_eulers_for_z_and_x_alignment(direction_vector_current_sampling_point, direction_vector)
    coordinates_and_direction_vectors_of_all_sampling_points = [[[round(point[0] + direction_vector_current_sampling_point[0],2),
                                                                  round(point[1] + direction_vector_current_sampling_point[1],2),
                                                                  round(point[2] + direction_vector_current_sampling_point[2],2)],
                                                                 eulers_current_vector]]
    for i in range(1, sampling_rate_circle, 1):
        direction_vector_current_sampling_point = rotate_vector_around_axis(direction_vector_current_sampling_point, normalized_direction_vector, sampling_radians)
        eulers_current_vector = all_eulers_for_z_and_x_alignment(direction_vector_current_sampling_point, direction_vector)
        coordinates_and_direction_vectors_of_all_sampling_points.append([[round(point[0] + direction_vector_current_sampling_point[0],2),
                                                                          round(point[1] + direction_vector_current_sampling_point[1],2),
                                                                          round(point[2] + direction_vector_current_sampling_point[2],2)],
                                                                         eulers_current_vector])
    return coordinates_and_direction_vectors_of_all_sampling_points


def sample_tube(contour, sampling_rate_contour, tube_radius, sampling_rate_circle):
    resampled_contour = contour_to_spline_interpolated_contour(contour, sampling_rate_contour)
    sampling_points_on_tube = []
    for point in range(0,len(resampled_contour)):
        if point == 0:
            sampling_points_on_tube += sample_tube_at_one_position(resampled_contour[0],
                                                                   geo_vector(resampled_contour[0], resampled_contour[1]),
                                                                   tube_radius,
                                                                   sampling_rate_circle)
        else:
            sampling_points_on_tube += sample_tube_at_one_position(resampled_contour[point],
                                                                   geo_vector(resampled_contour[point - 1],
                                                                              resampled_contour[point]),
                                                                   tube_radius, sampling_rate_circle)
    return sampling_points_on_tube


def sample_object(object, sampling_rate_contour, tube_radius, sampling_rate_circle):
    sampling_points_in_object = []
    for contour in object:
        sampling_points_in_object.append(sample_tube(contour, sampling_rate_contour, tube_radius, sampling_rate_circle))
    return sampling_points_in_object


def write_table_line(particle_tag,
                     tdrot,
                     tilt,
                     narot,
                     y_min_tilt,
                     y_max_tilt,
                     x_coord,
                     y_coord,
                     z_coord,
                     tube_number):
    #table_line = f"{particle_tag}\t0\t0\t0\t0\t0\t{tdrot}\t{tilt}\t{narot}\t0\t0\t0\t1\t{y_min_tilt}\t{y_max_tilt}\t0\t0\t0\t0\t0\t0\t0\t0\t{x_coord}\t{y_coord}\t{z_coord}\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t{tube_number}"
    table_line = f"{particle_tag} 1 0 0 0 0 {tdrot} {tilt} {narot} 0 0 0 1 {y_min_tilt} {y_max_tilt} 0 0 0 0 0 0 0 {tube_number} {x_coord} {y_coord} {z_coord} 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0"
    return table_line


def matlab_table_from_sampled_object(sampled_object, y_min_tilt, y_max_tilt, tube_number=1):
    particle_tag = 1
    matlab_table=""
    for sampled_tube in sampled_object:
        for sampling_point_on_tube in sampled_tube:
            coords = sampling_point_on_tube[0]
            eulers = sampling_point_on_tube[1]
            matlab_table += write_table_line(particle_tag, eulers[0], eulers[1], eulers[2], y_min_tilt, y_max_tilt, coords[0], coords[1], coords[2],tube_number) + "\n"
            particle_tag += 1
        tube_number += 1
    return matlab_table


def write_matlab_table_to_file(matlab_table, file_name):
    p = open(file_name, "w+")
    p.write(matlab_table)
    p.close()
    return []


def write_imod_txt_line(x_coord, y_coord, z_coord, tube_number):
    table_line = f"{tube_number}\t{x_coord}\t{y_coord}\t{z_coord}"
    return table_line


def imod_txt_from_sampled_object(sampled_object):
    imod_txt =  ""
    tube_number = 1
    for sampled_tube in sampled_object:
        for sampling_point_on_tube in sampled_tube:
            coords = sampling_point_on_tube[0]
            imod_txt += write_imod_txt_line(coords[0], coords[1], coords[2], tube_number) + "\n"
        tube_number += 1
    return imod_txt


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='samples the contours of an object assuming they represent filaments/tubes'
                                                 'output is in matlab .tbl format')
    parser.add_argument('object',
                        type=str,
                        help='object to be sampled in imod model2point -c .txt format')
    parser.add_argument('--tube_length_sampling_point_distance',
                        default=1,
                        type=float,
                        help='defines defines the distance between sampling points along the filament/tube axis'
                             'default = 1')
    parser.add_argument('--tube_radius',
                        default=0.1,
                        type=float,
                        help='defines in which radius (in pixels) around the central contour sampling will be performed.'
                             'default assumes simple filament')
    parser.add_argument('--tube_circumference_sampling_points',
                        default=1,
                        type=int,
                        help='defines how many regularly spaced points on an assumed tube circumference'
                             'are sampled. '
                             'default = 1')
    parser.add_argument('--min_y_tilt',
                        default=-60,
                        type=float,
                        help='defines the highest tilt angle on the minus side for the dynamo table.'
                             'default = -60')
    parser.add_argument('--max_y_tilt',
                        default=60,
                        type=float,
                        help='defines the highest tilt angle on the plus side for the dynamo table.'
                             'default = 60')
    parser.add_argument('--starting_tube_number',
                        default=1,
                        type=int,
                        help='defines the tube number (column 43) of the first contour in the object.'
                             'default = 1')
    parser.add_argument('--output',
                        default='output.tbl',
                        type=str,
                        help='defines the name of the .tbl that is written out.'
                             'default = output.tbl')

    args = parser.parse_args()
    write_matlab_table_to_file(matlab_table_from_sampled_object(sample_object(args.object,
                                                                              args.tube_length_sampling_point_distance,
                                                                              args.tube_radius,args.tube_circumference_sampling_points),
                                                                args.min_y_tilt,
                                                                args.max_y_tilt,
                                                                args.starting_tube_number),
                               args.output)

    imod_txt_from_sampled_object(sample_object(args.object,
                                               args.tube_length_sampling_point_distance,
                                               args.tube_radius,args.tube_circumference_sampling_points))







    #point1 = [0, 0, 0]
    #point2 = [2, 0, 0]
    #point3 = [0, 2, 0]
    #point4 = [0, 0, 2]
    #point5 = [2, 2, 2]

    #contour1 = [point1, point2]
    #contour2 = [point1, point3]
    #contour3 = [point1, point4]
    #contour4 = [point1, point5]

    #object1 = [contour1, contour2]



    #print(write_table_line(1, 10, 20, 30, -60, 60, 1, 2, 3, 5))
    #object1 = file_to_model("fils1Cleaned.txt")
    #print(matlab_table_from_sampled_object(sampled_object1, -60, 60, 5))
    #write_matlab_table_to_file(matlab_table_from_sampled_object(sample_object(object1,1,7,20), -60, 60, 5), "test_table.txt")



    #print(contour_to_spline_interpolated_contour(contour1,0.2))

    #print(sample_tube_at_one_position(point1, point2,2,2))
    #print(sample_tube(contour1,1,2,2))
    #print(sample_tube(contour2,1,2,2))
    #print(sample_object(object1,1,2,2))
    #print(sample_tube(contour3,1,2,2))
    #print(sample_tube(contour4,1,2,1))


    #erster euler = Winkel zwischen Projection der Contour auf die xy-Ebene und der x-Achse
    #zweiter_euler = Winkel zwischen Projection der Contour auf die xy-Ebene und Richtungsvektor des Samplingpoints
