cdef class VTFInputFile:
    cdef list configs

    cpdef list get_config_positions(self, int step)
