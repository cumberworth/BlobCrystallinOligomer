cdef extern from "BlobCrystallinOligomer/particle.h" namespace "particle":
    cdef cppclass Particle:
        int get_index();
        int get_type();
