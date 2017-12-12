from config_c cimport Config as Config_c
from random_gens_c cimport RandomGens as RandomGens_c

cdef class Config:
    cdef Config_c* config_c
    cdef RandomGens_c* random_gens_c
    cdef int num_monomers
    cpdef calc_dist_pairs(self, int particle_i)
