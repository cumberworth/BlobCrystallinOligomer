from libcpp.vector cimport vector

from particle_c cimport Particle
from shared_types_c cimport *

cdef extern from "BlobCrystallinOligomer/monomer.h" namespace "monomer":
    cdef cppclass Monomer:
        int get_index();
        int get_num_particles();
        Particle& get_particle(int particle_i);
        distT get_radius();
