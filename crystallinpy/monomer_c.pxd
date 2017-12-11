from libcpp.vector cimport vector

from shared_types_c cimport *

cdef extern from "../include/BlobCrystallinOligomer/monomer.h" namespace "monomer":
    cdef cppclass Monomer:
        int get_index();
        int get_num_particles();
        distT get_radius();
