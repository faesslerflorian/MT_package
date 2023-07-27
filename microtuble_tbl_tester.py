#import microtubule_tbl_from_imod_object as mtl
#import microtubule_tbl_from_imod_object_experimental as mtl
#import microtubule_tbl_from_imod_object_ex2 as mtl
import microtubule_tbl_from_imod_object_Z_align_with_center as mtl

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
#object1 = mtl.file_to_model("mapA5_area1__ts_005.txt")
#mtl.write_matlab_table_to_file(mtl.matlab_table_from_sampled_object(mtl.sample_object(object1,1,0.01,4), -65.8, 56.5, 1), "1_px_spaced_mapA5_area1__ts_005_center.tbl")
#object2 = mtl.file_to_model("mapA6_area1__ts_002.txt")
#mtl.write_matlab_table_to_file(mtl.matlab_table_from_sampled_object(mtl.sample_object(object2,1,0.01,4), -64.1, 55.2, 4), "1_px_spaced_mapA6_area1__ts_002_center.tbl")
#object3 = mtl.file_to_model("mapA11_area1__ts_003.txt")
#mtl.write_matlab_table_to_file(mtl.matlab_table_from_sampled_object(mtl.sample_object(object3,1,0.01,4), -67.3, 54.2, 43), "1_px_spaced_mapA11_area1__ts_003_center.tbl")
object4 = mtl.file_to_model("TS_14_intercisternal_filaments_bin8.txt")
mtl.write_matlab_table_to_file(mtl.matlab_table_from_sampled_object(mtl.sample_object(object4,2,0.01,1), -42, 56.1, 1), "TS_14_intercisternal_filaments_bin8.tbl")