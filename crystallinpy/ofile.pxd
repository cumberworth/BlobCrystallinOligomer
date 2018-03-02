cdef class VTFOutputFile:
    cdef object filename
    cdef object header
    cdef object ofile

    cpdef void write_config_positions(self, list positions)
