cdef class VTFInputFile:
    def __init__(self, filename):
        with open(filename) as inp:
            lines = inp.readlines()

        # Extract positions for each step to list of list of lists
        parsing_header = True
        self.configs = []
        for line in lines:
            words = line.split()
            if len(words) == 0:
                continue
            elif words[0] == 't':
                self.configs.append([])
                if parsing_header:
                    parsing_header = False
            else:
                if parsing_header:
                    continue
                else:
                    pos = [float(a) for a in words]
                    self.configs[-1].append(pos)

    @property
    def num_configs(self):
        return len(self.configs)

    cpdef list get_config_positions(self, int step):
        return self.configs[step]
