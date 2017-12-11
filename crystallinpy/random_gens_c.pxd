cdef extern from "../include/BlobCrystallinOligomer/random_gens.h" namespace "random_gens":
    cdef cppclass RandomGens:
        RandomGens() except +
