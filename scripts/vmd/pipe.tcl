# Shitty pipe

set libdir [lindex $argv 0]
set filebase [lindex $argv 1]

set num_vecs 3

source $libdir/libcrystallin.tcl

set system [mol new $filebase.vtf]
set filebase ${filebase}_pipe
mol delrep 0 0
create_reps
animate read vtf $filebase.vtf waitfor all $system
animate read vtf $filebase.vtf waitfor all $system
axes location off
display projection orthographic
mol top $system

trace variable vmd_frame(0) w update_frame_trace

animate speed 0.1
animate forward
