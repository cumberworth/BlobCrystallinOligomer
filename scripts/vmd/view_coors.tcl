# Requires the variables filebase and libdir have been set

source $libdir/libcrystallin.tcl

set num_vecs 2

set vecs_raw [load_matrix_as_lists $filebase.patch]
set vecs [unpack_vecs $vecs_raw $num_vecs]

set system [mol new $filebase.vtf waitfor all]
mol delrep 0 $system
create_reps
#mol addfile $filebase.vtf type vtf waitfor all
axes location off
display projection orthographic
mol top $system
trace variable vmd_frame(0) w draw_patch_vectors_trace

animate goto start
