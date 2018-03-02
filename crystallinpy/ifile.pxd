cdef class VTFInputFile:
    cdef list configs
    cdef object header

    cpdef list get_config_positions(self, int step)
