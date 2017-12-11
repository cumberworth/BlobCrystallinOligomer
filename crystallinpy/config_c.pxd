from libcpp.vector cimport vector

from ifile_c cimport MonomerData
from monomer_c cimport Monomer
from random_gens_c cimport RandomGens
from shared_types_c cimport *

cdef extern from "../include/BlobCrystallinOligomer/config.h" namespace "config":
    cdef cppclass Config:
        Config(
                vector[MonomerData] monomers,
                RandomGens& random_number,
                distT box_len,
                distT radius) except +
        Monomer& get_monomer(int monomer_index)
        int get_num_particles()
        distT get_box_len()
        distT get_radius()
        update_config_positions(vector[vector[double]])
