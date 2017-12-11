# distutils: language = c++
# distutils: sources = ../src/config.cpp

from libcpp.string cimport string
from libcpp.vector cimport vector

from config_c cimport Config as Config_c
from ifile_c cimport InputConfigFile as InputConfigFile_c
from ifile_c cimport MonomerData as MonomerData_c
from monomer_c cimport Monomer as Monomer_c
from random_gens_c cimport RandomGens as RandomGens_c
from shared_types_c cimport *

from monomer cimport Monomer
from monomer import Monomer

cdef class Config:
    cdef Config_c* config_c
    cdef RandomGens_c* random_gens_c

    def __init__(self, filename):
        pass

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

    def update_config_positions(positions):
        self.config_c.update_config_positions(positions)
