# Functions for visualizing patchy particle simulations
#
# system must be set to the relevant mol id

# Color scheme
set colors(0) 10
set colors(1) 11
set colors(2) 12
set colors(acd) 23
set colors(ntd) 31
set colors(blob) 1

proc load_matrix_as_lists {filename} {
    # Load a matrix as a list of lists
    set f [open $filename r]
    set raw [read $f]
    close $f
    set lines [split $raw "\n"]
    set last_index [expr [llength $lines] - 1]
    set lines [lreplace $lines $last_index $last_index]

    return $lines
}

proc unpack_vecs {raw_mat num_vecs} {
    # Take a list of lists matrix and package the vectors into a further list
    set vecs {}
    foreach step $raw_mat {
        set step_vecs {}
        for {set pi 0} {$pi != [llength $step]} {incr pi [expr 3*$num_vecs]} {
            set part {}
            set j $pi
            for {set vi 0} {$vi != $num_vecs} {incr vi} {
                set vec {}
                for {set k 0} {$k != 3} {incr k} {
                    lappend vec [lindex $step $j]
                    incr j
                }
                lappend part $vec
            }
            lappend step_vecs $part
        }
        lappend vecs $step_vecs
    }

    return $vecs
}

proc create_reps {} {
    global system
    global colors
    mol rep vdw
    mol addrep $system
    mol modselect 0 $system "type 0 1"
    mol modcolor 0 $system ColorID $colors(acd)

    mol addrep $system
    mol modselect 1 $system "type 2 3"
    mol modcolor 1 $system ColorID $colors(ntd)

    mol addrep $system
    mol modselect 2 $system "type 4"
    mol modcolor 2 $system ColorID $colors(blob)
}

# Not sure how to pass arguments properly with callbacks, so just use globals
# It's what they do in the examples in the VMD docs

proc draw_3d_vector {origin vector color} {
    # Draw vector from origin
    global system
    graphics $system color $color
    set end [vecadd $origin $vector]
    set middle [vecadd $origin [vecscale 0.8 [vecsub $end $origin]]]
    graphics $system cylinder $origin $middle radius 1.0 resolution 10
    graphics $system cone $middle $end radius 3.0 resolution 10
}

proc draw_patch_vectors {vecs num_vecs} {
    # Draw patch vectors for current frame
    # Must clear previous first
    global system
    global colors
    set frame [molinfo $system get frame]
    set num_particles [molinfo $system get numatoms]
    for {set i 0} {$i != $num_particles} {incr i} {
        set d [atomselect $system "index $i" frame $frame]
        set d_coors [lindex [$d get {x y z}] 0]
        for {set j 0} {$j != $num_vecs} {incr j} {
            set d_vec [lindex [lindex [lindex $vecs $frame] $i] $j]
            set vector [vecscale $d_vec 10.8]
            draw_3d_vector $d_coors $vector $colors($j)
        }
    }
}

proc draw_patch_vectors_trace {args} {
    global system
    global vecs
    global num_vecs
    graphics $system delete all

    # Redraw pbc box
    pbc box_draw -center origin
    draw_patch_vectors $vecs $num_vecs
}

proc update_frame {} {
    # Load new configuration and delete previous
    global system
    global filebase
    global num_vecs
    animate delete beg 0 $system

    # Save visulation state
    foreach mol [molinfo list] {
        set viewpoints($mol) [molinfo $mol get {
            center_matrix rotate_matrix scale_matrix global_matrix}]
    }
    animate read vtf $filebase.vtf waitfor all $system
    animate read vtf $filebase.vtf waitfor all $system

    # Return to previous visulation state
    foreach mol [molinfo list] {
        molinfo $mol set {center_matrix rotate_matrix scale_matrix
            global_matrix} $viewpoints($mol)
    }

    set vecs_raw [load_matrix_as_lists $filebase.patch]
    set vecs [unpack_vecs $vecs_raw $num_vecs]
    lappend vecs [lindex $vecs 0]
    graphics $system delete all

    # Redraw pbc box
    pbc box_draw -center origin
    draw_patch_vectors $vecs $num_vecs
}

proc update_frame_trace {args} {
    update_frame
}

proc align {} {
    global system
    global vecs
    set refsel [atomselect $system all frame 0]
    set a_vecs [list [lindex $vecs 0]]
    for {set i 1} {$i != [molinfo 0 get numframes]} {incr i} {
        set cursel [atomselect $system all frame $i]
        set alignmat [measure fit $cursel $refsel]
        $cursel move $alignmat

        set frame_vecs [lindex $vecs $i]
        set a_frame_vecs [list]
        foreach pvecs $frame_vecs {
            set a_pvecs [list]
            foreach vec $pvecs {
                lappend a_pvecs [vectrans $alignmat $vec]
            }
            lappend a_frame_vecs $a_pvecs
        }
        lappend a_vecs $a_frame_vecs
    }
    set vecs $a_vecs
    display resetview
}
