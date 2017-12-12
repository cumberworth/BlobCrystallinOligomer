cdef extern from "BlobCrystallinOligomer/shared_types.h" namespace "shared_types":
    ctypedef double distT;
    ctypedef double eneT;
    ctypedef unsigned long long int stepT;
    ctypedef double timeT;
    cdef enum CoorSet:
        current "shared_types::CoorSet::current",
        trial "shared_types::CoorSet::trial",
