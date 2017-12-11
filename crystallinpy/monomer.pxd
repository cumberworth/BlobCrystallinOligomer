from monomer_c cimport Monomer as Monomer_c

cdef class Monomer:
    cdef Monomer_c* monomer_c

    @staticmethod
    cdef Monomer create(Monomer_c* monomer_c)
