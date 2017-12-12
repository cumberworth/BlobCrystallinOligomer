from libcpp.string cimport string
from libcpp.vector cimport vector

from shared_types_c cimport *

cdef extern from "BlobCrystallinOligomer/ifile.h" namespace "ifile":
    cdef struct MonomerData:
        pass

    cdef cppclass InputConfigFile:
        InputConfigFile(string filename);
        vector[MonomerData] get_monomers();
        distT get_box_len();
        distT get_radius();
