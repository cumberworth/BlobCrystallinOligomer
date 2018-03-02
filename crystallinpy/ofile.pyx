cdef class VTFOutputFile:
    def __init__(self, filename, header):
        self.filename = filename
        self.header = header

        self.ofile = open(filename, 'w')
        self.ofile.write(header)

    cpdef void write_config_positions(self, list positions):
        self.ofile.write('t\n')
        for pos in positions:
            self.ofile.write('{} {} {}\n'.format(pos[0], pos[1], pos[2]))

        self.ofile.write('\n')
