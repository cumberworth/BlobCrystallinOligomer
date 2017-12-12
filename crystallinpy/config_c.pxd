from libcpp.vector cimport vector

from ifile_c cimport MonomerData
from monomer_c cimport Monomer
from particle_c cimport Particle
from random_gens_c cimport RandomGens
from shared_types_c cimport *

cdef extern from "BlobCrystallinOligomer/config.h" namespace "config":
    cdef cppclass Config:
        Config(
                vector[MonomerData] monomers,
                RandomGens& random_number,
                distT box_len,
                distT radius) except +
        Monomer& get_monomer(int monomer_index)
        int get_num_particles()
        int get_num_monomers()
        distT get_box_len()
        distT get_radius()
        void update_config_positions(vector[vector[double]] positions)
        distT calc_dist(
                Particle& particle1,
                CoorSet& coorset1,
                Particle& particle2,
                CoorSet& coorset2);
