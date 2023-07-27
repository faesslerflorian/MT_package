import microtubule_tbl_from_imod_object_Z_align_with_center as mtl

# mtl.sample_object(imported imod object
#                   number of sampling points in a radius around the central axis, 1 for your filaments
#                   radius which is sampled around the central axis, it can't be 0, stick to 0.01 for your purposes
#                   sampling distance in pixels (same as in imported imod object) along the axis, choose whatever makes sense in your case)
# mtl.matlab_table_from_sampled_object(mtl.sample_object ouput
#                                      highest tilt on the minus side, if you don't do STA in dynamo just put -60
#                                      highest tilt on the plus side, if you don't do STA in dynamo just put +60
#                                      starting number for filaments in this tomo, if you don't do STA in dynamo just put 1)
# mtl.write_matlab_table_to_file(mtl.matlab_table_from_sampled_object output
#                                desired name of the table file)


object1 = mtl.file_to_model("mapA5_area1__ts_005.txt")
mtl.write_matlab_table_to_file(mtl.matlab_table_from_sampled_object(mtl.sample_object(object1,1,0.01,4), -65.8, 56.5, 1), "1_px_spaced_mapA5_area1__ts_005_center.tbl")

object2 = mtl.file_to_model("mapA6_area1__ts_002.txt")
mtl.write_matlab_table_to_file(mtl.matlab_table_from_sampled_object(mtl.sample_object(object2,1,0.01,4), -64.1, 55.2, 4), "1_px_spaced_mapA6_area1__ts_002_center.tbl")

object3 = mtl.file_to_model("mapA11_area1__ts_003.txt")
mtl.write_matlab_table_to_file(mtl.matlab_table_from_sampled_object(mtl.sample_object(object3,1,0.01,4), -67.3, 54.2, 43), "1_px_spaced_mapA11_area1__ts_003_center.tbl")