from libcpp.string cimport string
from libcpp.vector cimport vector
cimport numpy as np

from config_c cimport Config as Config_c
from ifile_c cimport InputConfigFile as InputConfigFile_c
from ifile_c cimport MonomerData as MonomerData_c
from monomer_c cimport Monomer as Monomer_c
from particle_c cimport Particle as Particle_c
from random_gens_c cimport RandomGens as RandomGens_c
from shared_types_c cimport *

import numpy as np

from monomer cimport Monomer
from monomer import Monomer

DTYPE = np.float
ctypedef np.float_t DTYPE_t

cdef class Config:
    def __init__(self, filename):
        self.num_monomers = self.config_c.get_num_monomers()

    def __cinit__(self, filename):
        cdef string filename_c = filename.encode('UTF-8')
        cdef InputConfigFile_c* inputconfigfile = new InputConfigFile_c(filename_c)
        cdef vector[MonomerData_c] monomers = inputconfigfile.get_monomers()
        cdef distT box_len = inputconfigfile.get_box_len()
        cdef distT radius = inputconfigfile.get_radius()
        del inputconfigfile
        self.random_gens_c = new RandomGens_c()
        self.config_c = new Config_c(monomers, self.random_gens_c[0], box_len, radius)
 
    def __dealloc__(self):
        del self.config_c

    @property
    def num_monomers(self):
        return self.num_monomers

    def get_monomer(self, int monomer_i):
        cdef Monomer_c* monomer_c = &self.config_c.get_monomer(monomer_i)
        monomer = Monomer.create(monomer_c)
        return monomer

    def get_num_particles(self):
        return self.config_c.get_num_particles()

    def get_box_len(self):
        return self.config_c.get_box_len()

    def get_radius(self):
        return self.config_c.get_radius()

    def update_config_positions(self, positions):
        cdef vector[vector[double]] positions_v
        for i, pos in enumerate(positions):
            positions_v.push_back([])
            for x in pos:
                positions_v[i].push_back(x)

        self.config_c.update_config_positions(positions_v)

    cpdef calc_dist_pairs(self, int particle_i):
        cdef dists = np.zeros(self.num_monomers*(self.num_monomers - 1)/2, dtype=DTYPE)
        cdef int i, j, pair_i
        cdef Monomer_c* m1
        cdef Monomer_c* m2
        cdef Particle_c* p1
        cdef Particle_c* p2
        cdef distT dist
        cdef CoorSet current = CoorSet.current
        pair_i = 0
        for i in range(self.num_monomers):
            m1 = &self.config_c.get_monomer(i)
            p1 = &m1.get_particle(particle_i)
            for j in range(i + 1, self.num_monomers):
                m2 = &self.config_c.get_monomer(j)
                p2 = &m2.get_particle(particle_i)
                dist = self.config_c.calc_dist(
                        p1[0], current,
                        p2[0], current)
                dists[pair_i] = dist
                pair_i += 1

        return dists
