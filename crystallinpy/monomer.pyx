from monomer_c cimport Monomer as Monomer_c

cdef class Monomer:

    @staticmethod
    cdef Monomer create(Monomer_c* monomer_c):
        obj = <Monomer>Monomer.__new__(Monomer)
        obj.monomer_c = monomer_c
        return obj

    def get_index(self):
        return self.monomer_c.get_index()

    def get_num_particles(self):
        return self.monomer_c.get_num_particles()

    def get_radius(self):
        return self.monomer_c.get_radius()
